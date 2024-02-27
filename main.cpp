#include <iostream>
#include <iomanip>
#include "convert.hpp"

int main(int argc, char **argv)
{
    if (argc != 4)
    {
        std::cerr << "Usage: " << argv[0] << " <lat> <lon> <alt>" << std::endl;
        return -1;
    }

    double wLat = atof(argv[1]);
    double wLon = atof(argv[2]);
    double wAlt = atof(argv[3]);
    double sLat = 0.0, sLon = 0.0, sAlt = 0.0;
    double sX = 0.0, sY = 0.0;

    std::cout << "Coordinate in WGS84: (" << wLat << ", " << wLon << ", " << wAlt << ")\n";

    Converter converter;
    std::tie(sLat, sLon, sAlt) = converter.WGS84ToSK42(wLat, wLon, wAlt);
    std::tie(sX, sY) = converter.ProjectSK42(sLat, sLon);

    std::cout << std::setprecision(18) << "Coordinate in SK42: (" << sLat << ", " << sLon << ", " << sAlt << ")\n"
              << "Coordinate in SK42 Project: (" << sX << ", " << sY << ")" << std::endl;

    /*
     * Coordinate in WGS84: 113.95544, 34.13451, 100

     * From https://epsg.io/srs/transform/113.95544,34.13451,100.0.json?key=default&s_srs=4326&t_srs=4284
     * Coordinate in SK42: 113.95492242191457, 34.1343303935188, 150.24392636213452
     * Our coordinate:     113.95492242597204, 34.1343303934879, 150.24400700908154

     * From https://sk42.org/?direction=toSK42&convstr=34.13451%2C113.95544
     * Coordinate in SK42 project:     3782576            , 19772606
     * Our coordinate in SK42 project: 3782575.55743047688, 19772600.2419115938
     */

    return 0;
}