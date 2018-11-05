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
    std::cout << "\n";
}

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
        geometry.push_back(tool(20));   // TODO simple circle

        std::vector<std::vector<point_2>> stock;
        {
            stock.emplace_back();
            auto& path = stock.back();
            path.push_back({-30, -30});
            path.push_back({30, -30});
            path.push_back({30, 30});
            path.push_back({-30, 30});
        }

        double tool_diameter = 5.0;

        cl::Clipper clipper;
        clipper.AddPaths(scale_paths(stock), cl::ptSubject, true);
        clipper.AddPaths(scale_paths(geometry), cl::ptClip, true);
        cl::Paths material_to_remove;
        clipper.Execute(cl::ctDifference, material_to_remove);

        output_path(material_to_remove);

        cl::ClipperOffset co;
        co.AddPaths(material_to_remove, cl::jtRound, cl::etClosedPolygon);
        double step = tool_diameter / 2.0;// (180.0 / vm["engagement"].as<double>());
        double offset = -step;
        double target_engagement = vm["engagement"].as<double>();

        std::vector<cl::Paths> levels;
        while (true) {
            cl::Paths level;
            co.ArcTolerance = 0.1 * CLIPPER_SCALE;            
            co.Execute(level, offset * CLIPPER_SCALE);

            if (level.empty())
                break;

            levels.insert(begin(levels), level);

            offset -= step;
        }

        auto translate = [](const cl::Path& path, cl::IntPoint p) {
            cl::Path t;
            for (auto& px : path)
                t.emplace_back(px.X + p.X, px.Y + p.Y);
            return t;
        };
        auto lerp = [](const cl::IntPoint& P0, const cl::IntPoint& P1, double t) {
            auto p0 = unscale_point(P0);
            auto p1 = unscale_point(P1);
            return scale_point({ (1-t)*p0.x + t*p1.x, (1-t)*p0.y + t*p1.y });
        };
        auto distance = [](const point_2& a, const point_2& b) {
            return std::sqrt(std::pow(b.x - a.x, 2) + std::pow(b.y - a.y, 2));
        };
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

        auto tool_model = scale_path(tool(tool_diameter/2.0));
        // TODO contour offset material to remove!
        for (auto& level : levels) {
            for (auto& path : level) {
                auto pn1 = path.back();
                //std::cout << "M " << r6(first.x) << " " << r6(first.y) << " ";
                for (auto& pn : path) {
                    auto d = distance(unscale_point(pn1), unscale_point(pn));
                    auto steps = d / .1;
                    if (steps < 1) steps = 1;

                    auto step = 1.0/steps;
                    for (double t = 0; t <= 1; t += step) {
                        // line from pn-1 - pn
                        auto tp = translate(tool_model, lerp(pn1, pn, t));

                        cl::Clipper clipper;
                        clipper.AddPaths(material_to_remove, cl::ptClip, true);
                        clipper.AddPath(tp, cl::ptSubject, false);

                        cl::PolyTree pt;
                        clipper.Execute(cl::ctIntersection, pt);

                        auto l0 = path_length({tp});
                        cl::Paths open;
                        OpenPathsFromPolyTree(pt, open);

                        double engagement = -1;
                        if (!open.empty()) {
                            auto l1 = path_length(open);
                            engagement = 360.0 * (l1 / l0);
                            
                            if (std::abs(engagement - target_engagement) > 10.0) {
                                //auto dx = pn.x - pn1.x;
                                //auto dy = pn.y - pn1.y;
                                //auto p3 = vector_2(-dy, dx);

                            }
                        }
                        //std::cout << "L " << r6(pn.x) << " " << r6(pn.y) << " ";

                        {
                            cl::Clipper clipper;
                            clipper.AddPaths(material_to_remove, cl::ptSubject, true);
                            clipper.AddPath(tp, cl::ptClip, true);

                            material_to_remove.clear();
                            clipper.Execute(cl::ctDifference, material_to_remove);
                        }
                    }
                    pn1 = pn;
                }
            }
        }

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

