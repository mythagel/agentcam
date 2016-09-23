/*
 * Copyright (C) 2014  Nicholas Gill
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include <iostream>
#include <iomanip>
#include <string>
#include <stdexcept>
#include <boost/program_options.hpp>
#include <vector>
#include <algorithm>

#include "clipper.hpp"

namespace po = boost::program_options;
namespace cl = ClipperLib;

static const double CLIPPER_SCALE = 1e12;

void print_exception(const std::exception& e, int level = 0) {
    std::cerr << std::string(level, ' ') << "error: " << e.what() << '\n';
    try {
        std::rethrow_if_nested(e);
    } catch(const std::exception& e) {
        print_exception(e, level+1);
    } catch(...) {
    }
}

inline std::string r6(double v) {
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(6) << v;
    auto s = ss.str();
    
    s.erase(s.find_last_not_of('0') + 1, std::string::npos);
    if(s.back() == '.') s.pop_back();
    return s;
}


/* staticcam
 * 3d model -> 2d slices at DOC
 * generate 2d slice
 * tool as 2d projection
 * genetic algorithm minimise machining time
 *
 * tool is a circle at the effective diameter
 *
 * fixed paramters
 * spindle speed range
 * spindle torque
 * tool stickout
 * tool diameter
 * machine rigidity
 * */

const double PI = 3.14159265359;

struct point_2 {
    double x;
    double y;
    point_2()
     : x(), y() {
    }
    point_2(double x, double y)
     : x(x), y(y) {
    }
};

std::vector<point_2> tool(double r, unsigned segments = 64) {
    std::vector<point_2> polygon;
    polygon.reserve(segments);

    double delta_theta = (2*PI) / segments;
    double theta = 0;
    for (unsigned i = 0; i < segments; ++i, theta += delta_theta)
        polygon.emplace_back(std::cos(theta)*r, std::sin(theta)*r);

    return polygon;
}

/* staticcam
 * generate toolpath based on (approximately) static tool engagement angle.
 * 2d only input (atm)
 * prefer climb over conventional
 *
 * from current point, need to determine next point to move
 * have theta representing direction and dist representing how far to travel
 * input to this function is the stock and part 2d polygons
 *
 */

point_2 unscale_point(const cl::IntPoint& p) {
    return point_2(static_cast<double>(p.X) / CLIPPER_SCALE, static_cast<double>(p.Y) / CLIPPER_SCALE);
}
cl::IntPoint scale_point(const point_2& p) {
    return cl::IntPoint(p.x * CLIPPER_SCALE, p.y * CLIPPER_SCALE);
}
cl::Path scale_path(const std::vector<point_2>& path) {
    cl::Path scaled;
    scaled.reserve(path.size());
    for (auto& p : path)
        scaled.push_back(scale_point(p));
    return scaled;
}
cl::Paths scale_paths(const std::vector<std::vector<point_2>>& paths) {
    cl::Paths scaled;
    scaled.reserve(paths.size());
    for (auto& path : paths)
        scaled.push_back(scale_path(path));
    return scaled;
}

void output_path(const cl::Paths& paths) {
    for(auto& path : paths) {
        if (path.empty()) continue;
        auto first = unscale_point(path.front());
        std::cout << "M " << r6(first.x) << " " << r6(first.y) << " ";
        for(auto& point : path) {
            auto p = unscale_point(point);
            std::cout << "L " << r6(p.x) << " " << r6(p.y) << " ";
        }
        std::cout << "L " << r6(first.x) << " " << r6(first.y) << " ";
    }
}

struct state {
    point_2 current;
    double current_theta = -1;

    double tool_radius;
    double target_engagement;
    cl::Paths stock;
    cl::Paths model;


    struct candidate_t {
        double theta = 0;
        double r = 0;

        // Calculated
        cl::Paths tool;
        double engagement = -1;

        point_2 generate_point(const point_2& current) const {
            point_2 center;
            center.x = current.x + (std::cos(theta) * r);
            center.y = current.y + (std::sin(theta) * r);
            return center;
        }

        void generate_tool(const point_2& current, double tool_radius) {
            tool.clear();
            tool.emplace_back();
            auto& path = tool.back();

            const unsigned segments = 128;
            path.reserve(segments);
            
            auto center = generate_point(current);

            double delta_theta = (2.0*PI) / segments;
            double theta = 0;
            for (unsigned i = 0; i < segments; ++i, theta += delta_theta) {
                point_2 p;
                p.x = center.x + (std::cos(theta)*tool_radius);
                p.y = center.y + (std::sin(theta)*tool_radius);
                path.push_back(scale_point(p));
            }
        }

        void generate_engagement(const point_2& current, double tool_radius, const cl::Paths& stock) {
            generate_tool(current, tool_radius);

            cl::Clipper clipper;
            clipper.AddPaths(stock, cl::ptClip, true);
            clipper.AddPaths(tool, cl::ptSubject, false);

            cl::PolyTree pt;
            clipper.Execute(cl::ctIntersection, pt);

            auto path_length = [](const cl::Paths& paths) {
                double length = 0;
                if (paths.empty()) return length;
                for (auto& path : paths) {
                    if (path.empty()) return length;
                    auto p0 = unscale_point(path.front());
                    for (auto& point : path) {
                        auto p1 = unscale_point(point);
                        length += std::sqrt(std::pow(p1.x - p0.x, 2) + std::pow(p1.y - p0.y, 2));
                        p0 = p1;
                    }
                }
                return length;
            };

            auto l0 = path_length(tool);
            cl::Paths open;
            OpenPathsFromPolyTree(pt, open);
            if (!open.empty()) {
                auto l1 = path_length(open);
                engagement = 360.0 * (l1 / l0);
            } else {
                engagement = -1;
            }
        }
    };

    bool next(point_2& next_point) {
        /* generate candidate rays (theta / radius)
         * - potential next point is ray end point
         * - ugly brute force approach
         *
         * calculate tool engagement for each ray (end point)
         * if any tool point intersects the model, reduce radius until
         * it does not (or remove potential point if not possible)
         *
         * increase radius while cutter engagement remains (approximately) constant
         *
         * classify potential cut as climb, conventional, or slot
         * - reject conventional, slot sorted last
         *
         * sort closeness to current vector direction (prefer smooth movement)
         */

        // create candidate rays
        std::vector<candidate_t> candidates;
        unsigned rays = 16;
        double delta_theta = (2*PI) / rays;
        double theta = 0;
        for (unsigned ray = 0; ray < rays; ++ray) {
            candidate_t ci;
            ci.theta = theta;
            ci.r = 0.3;
            candidates.push_back(ci);
            theta += delta_theta;
        }

        // adjust / prune for model intersection
        for (auto& candidate : candidates) {

            // if tool intersects model, adjust by reducing radius if possible else remove from set
            auto tool_intersects_model = [&]() {
                candidate.generate_tool(current, tool_radius);

                cl::Clipper clipper;
                clipper.AddPaths(model, cl::ptSubject, true);
                clipper.AddPaths(candidate.tool, cl::ptClip, true);

                cl::Paths intersection;
                clipper.Execute(cl::ctIntersection, intersection);
                return !intersection.empty();
            };

            if (!tool_intersects_model()) {
                candidate.generate_engagement(current, tool_radius, stock);
                auto engagement = candidate.engagement;
                auto equiv = [&]{
                    return std::abs(engagement - candidate.engagement) < 0.1;
                };

                if (false) {
                    while (engagement > 0 && equiv()) {
                        std::cerr << "engagement: " << engagement << " radius: " << candidate.r << "\n";
                        candidate.r += 0.03;
                        candidate.generate_engagement(current, tool_radius, stock);
                    }
                    candidate.r -= 0.03;
                    // expand while engagement remains approx constant
                }
            }

            while (candidate.r > 0 && tool_intersects_model()) {
                candidate.r -= 0.06;
            }
        }

        candidates.erase(std::remove_if(begin(candidates), end(candidates), [](const candidate_t& t) { return t.r <= 0; }), end(candidates));

        // Generate cutter engagements
        for (auto& candidate : candidates) {
            candidate.generate_engagement(current, tool_radius, stock);
        }

        candidates.erase(std::remove_if(begin(candidates), end(candidates), [](const candidate_t& t) { return t.engagement <= 0; }), end(candidates));

        std::sort(begin(candidates), end(candidates), [this](const candidate_t& c0, const candidate_t& c1) {
            auto cmp = [this](const candidate_t& c) {
                auto dt = std::abs(c.theta - current_theta);
                return std::make_tuple(dt < PI/4, std::abs(c.engagement - target_engagement));
            };
            return cmp(c0) < cmp(c1);
        });

        for (auto& candidate : candidates) {
            std::cerr << "t: " << candidate.theta << " engagement: " << candidate.engagement << " target: " << target_engagement << "\n";
        }

        if (!candidates.empty()) {
            auto& candidate = candidates.front();
            next_point = candidate.generate_point(current);
            std::cerr << "T: " << candidate.theta << " ENGAGEMENT: " << candidate.engagement << " TARGET: " << target_engagement << "\n";

            cl::Clipper clipper;
            clipper.AddPaths(stock, cl::ptSubject, true);
            clipper.AddPaths(candidate.tool, cl::ptClip, true);

            stock.clear();
            clipper.Execute(cl::ctDifference, stock);

            //output_path(candidate.tool);

            current = next_point;
            current_theta = candidate.theta;
            return true;
        }
        return false;
    }
};

int main(int argc, char* argv[]) {
    po::options_description options("staticcam");
    std::vector<std::string> args(argv, argv + argc);
    args.erase(begin(args));

    options.add_options()
        ("help,h", "display this help and exit")
        ("engagement,e", po::value<double>()->default_value(120), "Cutter engagement angle (degrees)")
    ;

    try {
        po::variables_map vm;
        store(po::command_line_parser(args).options(options).run(), vm);

        if(vm.count("help")) {
            std::cout << options << "\n";
            return 0;
        }
        notify(vm);

        std::vector<std::vector<point_2>> geometry;
        geometry.push_back(tool(12.5));   // TODO simple circle

        std::vector<std::vector<point_2>> stock;
        {
            stock.emplace_back();
            auto& path = stock.back();
            path.push_back({-15, -15});
            path.push_back({15, -15});
            path.push_back({15, 15});
            path.push_back({-15, 15});
        }

        double tool_diameter = 5.0;

        state s;
        s.current = point_2(-15, -15);
        s.tool_radius = tool_diameter/ 2.0;
        s.target_engagement = vm["engagement"].as<double>();
        s.stock = scale_paths(stock);
        s.model = scale_paths(geometry);

        output_path(scale_paths(geometry));
        output_path(s.stock);

        point_2 point = s.current;
        std::cout << "M " << r6(point.x) << " " << r6(point.y) << " ";
        while (s.next(point)) {
            std::cout << "L " << r6(point.x) << " " << r6(point.y) << " ";
            std::cout << std::endl;
        }



        // Temp simple pocket

        /* agent must output path (list of points)
         * agent calculates next action, move to point!
         * how does agent decide how to move (i.e. which parameters are to be optimised?)
         * */

        // TODO
    } catch(const po::error& e) {
        print_exception(e);
        std::cout << options << "\n";
        return 1;
    } catch(const std::exception& e) {
        print_exception(e);
        return 1;
    }

    return 0;
}

