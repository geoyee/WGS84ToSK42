/*
 * Compute coordinate system transformation coefficients
 * For transformation from WGS-84 to Pulkovo 1942 (SK-42) and projection
 */

#ifndef CONVERT_HPP
#define CONVERT_HPP

#include <iostream>
#include <stdexcept>
#include <tuple>
#include <exception>
#include <cmath>

class Converter
{
public:
    // Constants
    const double M_PI = 3.141592653589794;
    const double M_Deg = M_PI / 180.0; // 0.0174532925199433
    const double M_Rad = 180.0 / M_PI; // 57.29577951308231
    const double M_EPS = 1e-12;
    /* From https://standartgost.ru/g/ГОСТ_Р_51794-2008 */
    // WGS84
    const double M_WGS84_a = 6378137.0;
    const double M_WGS84_f = 1.0 / 298.257223563;                                           // 0.0033528106647474805
    const double M_WGS84_b = M_WGS84_a * (1.0 - M_WGS84_f);                                 // 6356752.314245179
    const double M_WGS84_e = sqrt(1.0 - (M_WGS84_b * M_WGS84_b) / (M_WGS84_a * M_WGS84_a)); // 0.08181919084262149
    const double M_WGS84_e2 = M_WGS84_e * M_WGS84_e;                                        // 0.006694379990141316
    // Krassowsky
    const double M_Krassowsky_a = 6378245.0;
    const double M_Krassowsky_f = 1.0 / 298.3;                             // 0.003352329869259135
    const double M_Krassowsky_b = M_Krassowsky_a * (1.0 - M_Krassowsky_f); // 6356863.018773047
    const double M_Krassowsky_e =
        sqrt(1.0 - (M_Krassowsky_b * M_Krassowsky_b) / (M_Krassowsky_a * M_Krassowsky_a)); // 0.08181333401693065
    const double M_Krassowsky_e2 = M_Krassowsky_e * M_Krassowsky_e;                        // 0.006693421622965863

    explicit Converter(bool usingEPSG = true) : m_usingEPSG(usingEPSG) { }

    ~Converter() = default;

    std::tuple<double, double, double> WGS84ToSK42(double lat, double lon, double alt)
    {
        if (!checkLatLon(lat, lon))
        {
            std::invalid_argument("Lat or lon is out of range [-90, 90] or [-180, 180].");
        }
        double x = 0.0, y = 0.0, z = 0.0;
        std::tie(x, y, z) = BLHToXYZInWGS84(lat, lon, alt);
        std::tie(x, y, z) = m_usingEPSG ? WGS84ToSK42UsingEPSG(x, y, z) : WGS84ToSK42UsingGOST(x, y, z);
        double skLat = 0.0, skLon = 0.0, skAlt = 0.0;
        std::tie(skLat, skLon, skAlt) = XYZToBLHInKrassowsky(x, y, z);
        return std::make_tuple(skLat, skLon, skAlt);
    }

    std::tuple<double, double> ProjectSK42(double lat, double lon)
    {
        if (!checkLatLon(lat, lon))
        {
            std::invalid_argument("Lat or lon is out of range [-90, 90] or [-180, 180].");
        }
        return positiveGaussianInSK42(lat, lon);
    }

private:
    bool m_usingEPSG = true;

    inline double deg2Rad(double deg) noexcept
    {
        return deg * M_Deg;
    }

    inline double rad2Deg(double rad) noexcept
    {
        return rad * M_Rad;
    }

    inline bool checkLatLon(double lat, double lon) noexcept
    {
        return (lat >= -90.0 && lat <= 90.0) && (lon >= -180.0 && lon <= 180.0);
    }

    std::tuple<double, double, double> BLHToXYZInWGS84(double b, double l, double h) noexcept
    {
        double sinB = sin(deg2Rad(b)), cosB = cos(deg2Rad(b));
        double N = M_WGS84_a / sqrt(1.0 - M_WGS84_e2 * sinB * sinB);
        double x = (N + h) * cosB * cos(deg2Rad(l));
        double y = (N + h) * cosB * sin(deg2Rad(l));
        double z = ((1.0 - M_WGS84_e2) * N + h) * sinB;
        return std::make_tuple(x, y, z);
    }

    std::tuple<double, double, double> XYZToBLHInKrassowsky(double x, double y, double z) noexcept
    {
        double rl = atan2(y, x);
        double xy_hypot = hypot(x, y);
        double rb0 = 0.0;
        double rb = atan(z / xy_hypot);
        double N = 0.0;
        while (abs(rb - rb0) > M_EPS)
        {
            rb0 = rb;
            N = M_Krassowsky_a / sqrt(1.0 - M_Krassowsky_e2 * sin(rb0) * sin(rb0));
            rb = atan((z + M_Krassowsky_e2 * N * sin(rb0)) / xy_hypot);
        }
        N = M_Krassowsky_a / sqrt(1.0 - M_Krassowsky_e2 * sin(rb) * sin(rb));
        double h = 0.0;
        if (abs(rb) < M_PI / 4.0)
        {
            double R = hypot(xy_hypot, z);
            double phi = atan(z / xy_hypot);
            h = R * cos(phi) / cos(rb) - N;
        }
        else
        {
            h = z / sin(rb) - N * (1.0 - M_Krassowsky_e2);
        }
        return std::make_tuple(rad2Deg(rb), rad2Deg(rl), h);
    }

    /* Implemented according to ГОСТ 51794-2008 - 5.4 */
    std::tuple<double, double> positiveGaussianInSK42(double lat, double lon) noexcept
    {
        double bK = lat, lK = lon, bKRad = bK * M_Deg;
        int n = (int)((6.0 + lK) / 6.0); // Six-degree banding
        double l = (lK - abs(3.0 + 6.0 * (n - 1))) / M_Rad;
        double sinBk = sin(bKRad), sinBk2 = sinBk * sinBk, sinBk4 = sinBk2 * sinBk2, sinBk6 = sinBk4 * sinBk2;
        double l2 = l * l;
        double sx =
            6367558.4968 * bKRad -
            (sin(2 * bKRad)) *
                (16002.8900 + 66.9607 * sinBk2 + 0.3515 * sinBk4 -
                 l2 * (1594561.25 + 5336.535 * sinBk2 + 26.790 * sinBk4 + 0.149 * sinBk6 +
                       l2 * (672483.4 - 811219.9 * sinBk2 + 5420.0 * sinBk4 - 10.6 * sinBk6 +
                             l2 * (278194.0 - 830174.0 * sinBk2 + 572434.0 * sinBk4 - 16010.0 * sinBk6 +
                                   l2 * (109500.0 - 574700.0 * sinBk2 + 863700.0 * sinBk4 - 398600.0 * sinBk6)))));
        double sy = (5.0 + 10.0 * n) * pow(10.0, 5) +
                    l * cos(bKRad) *
                        (6378245.0 + 21346.1415 * sinBk2 + 107.1590 * sinBk4 + 0.5977 * sinBk6 +
                         l2 * (1070204.16 - 2136826.66 * sinBk2 + 17.98 * sinBk4 - 11.99 * sinBk6) +
                         l2 * (270806.0 - 1523417.0 * sinBk2 + 1327645.0 * sinBk4 - 21701.0 * sinBk6 +
                               l2 * (79690.0 - 866190.0 * sinBk2 + 1730360.0 * sinBk4 - 945460.0 * sinBk6)));
        return std::make_tuple(sx, sy);
    }

    /*
     * From https://standartgost.ru/g/ГОСТ_Р_51794-2008

        |X|   |            1  +3.1998*10^-6  -1.6968*10^-6|   |X|   | +25|
        |Y| = |-3.1998*10^-6              1              0| · |Y| - |-141|
        |Z|   |+1.6968*10^-6              0              1|   |Z|   | -80|
       SK-42                                                 PZ-90

        |X|                      |            1  +0.9696*10^-6              0|   |X|   |-1.10|
        |Y| = (1 + 0.12*10^-6) * |-0.9696*10^-6              1              0| · |Y| - |-0.30|
        |Z|                      |            0              0              1|   |Z|   |-0.90|
       PZ-90                                                                   WGS-84
     */
    inline std::tuple<double, double, double> WGS84ToSK42UsingGOST(double x, double y, double z) noexcept
    {
        auto _WGS84ToPZ90 = [](double x, double y, double z) -> std::tuple<double, double, double>
        {
            double k = 1.0 + 0.12e-6;
            double px = k * (x + 0.9696e-6 * y + 0.0) + 1.1;
            double py = k * (-0.9696e-6 * x + y + 0.0) + 0.3;
            double pz = k * (0.0 + 0.0 + z) + 0.9;
            return std::make_tuple(px, py, pz);
        };

        auto _PZ90ToSK42 = [](std::tuple<double, double, double> xyz) -> std::tuple<double, double, double>
        {
            double x = 0.0, y = 0.0, z = 0.0;
            std::tie(x, y, z) = xyz;
            double sx = (x + 3.1998e-6 * y - 1.6968e-6 * z) - 25.0;
            double sy = (-3.1998e-6 * x + y + 0.0) + 141.0;
            double sz = (1.6968e-6 * x + 0.0 + z) + 80.0;
            return std::make_tuple(sx, sy, sz);
        };

        return _PZ90ToSK42(_WGS84ToPZ90(x, y, z));
    }

    /*
     * From https://epsg.io/4284
     * Revision date: 2008-09-24

        TOWGS84[25,-141,-78.5,0,-0.35,-0.736,0]]
     */
    inline std::tuple<double, double, double> WGS84ToSK42UsingEPSG(double x, double y, double z) noexcept
    {
        double rY = -0.35 / 3600.0 * M_Deg, rZ = -0.736 / 3600.0 * M_Deg;
        double sx = (x - rZ * y + rY * z) - 25.0;
        double sy = (rZ * x + y - 0.0) + 141.0;
        double sz = (-rY * x + 0.0 + z) + 78.5;
        return std::make_tuple(sx, sy, sz);
    }
};

#endif // CONVERT_HPP