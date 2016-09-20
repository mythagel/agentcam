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
#include <boost/filesystem.hpp>
#include <vector>

extern "C" {
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"
}

#include "clipper.hpp"
#include "visilibity.hpp"

namespace po = boost::program_options;
namespace fs = boost::filesystem;
namespace cl = ClipperLib;
namespace vl = VisiLibity;

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

/* agentcam
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


struct path_scaling {
    const double scale = 1e12;

    cl::IntPoint scale_point(const point_2& p) const {
        return cl::IntPoint(p.x * scale, p.y * scale);
    }
    cl::Path scale_path(const std::vector<point_2>& path) const {
        cl::Path scaled;
        scaled.reserve(path.size());
        for (auto& p : path)
            scaled.push_back(scale_point(p));
        return scaled;
    }
    cl::Paths scale_paths(const std::vector<std::vector<point_2>>& paths) const {
        cl::Paths scaled;
        scaled.reserve(paths.size());
        for (auto& path : paths)
            scaled.push_back(scale_path(path));
        return scaled;
    }

    void output_path(const cl::Paths& paths) {
        auto unscale = [&](const cl::IntPoint& p) -> point_2 {
            return {static_cast<double>(p.X) / scale, static_cast<double>(p.Y) / scale};
        };
        for(auto& path : paths) {
            auto first = unscale(*path.begin());
            std::cout << "M " << r6(first.x) << " " << r6(first.y) << " ";
            for(auto& point : path) {
                auto p = unscale(point);
                std::cout << "L " << r6(p.x) << " " << r6(p.y) << " ";
            }
            std::cout << "L " << r6(first.x) << " " << r6(first.y) << " ";
        }
    }
};

vl::Polygon to_polygon(const std::vector<point_2>& points) {
    std::vector<vl::Point> polygon;
    polygon.reserve(points.size());
    for (auto& p : points)
        polygon.emplace_back(p.x, p.y);
    return { polygon };
}
vl::Environment to_environment() {
    return {};
}

class agent : private path_scaling {
private:
    cl::Path tool_p;
    cl::Paths geometry;
    point_2 current;

    void cut_path(const cl::IntPoint& p0, const cl::IntPoint& p1) {
        cl::Paths toolpath;
        MinkowskiSum(tool_p, {p0, p1}, toolpath, false);

        cl::Clipper clipper;
        clipper.AddPaths(geometry, cl::ptSubject, true);
        clipper.AddPaths(toolpath, cl::ptClip, true);

        geometry.clear();
        clipper.Execute(cl::ctDifference, geometry);
    }
public:
    agent (double tool_diameter, const std::vector<std::vector<point_2>>& geometry)
     : geometry(scale_paths(geometry)) {
        tool_p = scale_path(tool(tool_diameter/2.0));
    }
    point_2 step() {
        point_2 next;
        next = point_2{25, 25};
        cut_path(scale_point(current), scale_point(next));
        output_path(geometry);
        return next;
    }
};

int main(int argc, char* argv[]) {
    po::options_description options("agentcam");
    std::vector<std::string> args(argv, argv + argc);
    args.erase(begin(args));

    options.add_options()
        ("help,h", "display this help and exit")
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
        geometry.push_back(tool(25));   // TODO simple circle

        double tool_diameter = 4.0;
        auto bob = agent(tool_diameter, geometry);
        bob.step();

        // Temp simple pocket

        /* agent must output path (list of points)
         * agent calculates next action, move to point!
         * how does agent decide how to move (i.e. which parameters are to be optimised?)
         * */

        /* Part input as .off model, stock as .off model
         * stock - model == material to remove
         *
         * Make cuboid, height as depth of cut, width and depth greater than part dimensions
         * take intersection of 3d model using cork library, project as 2d plane
         * 2d plane slices represent material to remove per slice.
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

