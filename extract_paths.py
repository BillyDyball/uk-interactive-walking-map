#!/usr/bin/env python3
"""
Extract LDWA path data from search_by_path.php into GeoJSON.
Coordinates are converted from British National Grid (OSGB36) to WGS84 (lat/lng).
"""

import re
import json
import math
import html


def osgb36_to_wgs84(easting, northing):
    """Convert OSGB36 (British National Grid) to WGS84 lat/lng.
    Uses a simplified Helmert transformation — accurate to ~5m, fine for map display.
    """
    # Airy 1830 ellipsoid
    a, b = 6377563.396, 6356256.909
    F0 = 0.9996012717
    lat0 = math.radians(49)
    lon0 = math.radians(-2)
    N0, E0 = -100000, 400000
    e2 = 1 - (b * b) / (a * a)
    n = (a - b) / (a + b)

    lat = (northing - N0) / (a * F0) + lat0
    for _ in range(10):
        M = b * F0 * (
            (1 + n + 5 / 4 * n**2 + 5 / 4 * n**3) * (lat - lat0)
            - (3 * n + 3 * n**2 + 21 / 8 * n**3) * math.sin(lat - lat0) * math.cos(lat + lat0)
            + (15 / 8 * n**2 + 15 / 8 * n**3) * math.sin(2 * (lat - lat0)) * math.cos(2 * (lat + lat0))
            - (35 / 24 * n**3) * math.sin(3 * (lat - lat0)) * math.cos(3 * (lat + lat0))
        )
        if abs(northing - N0 - M) < 0.00001:
            break
        lat += (northing - N0 - M) / (a * F0)

    sin_lat = math.sin(lat)
    cos_lat = math.cos(lat)
    tan_lat = math.tan(lat)
    nu = a * F0 / math.sqrt(1 - e2 * sin_lat**2)
    rho = a * F0 * (1 - e2) / (1 - e2 * sin_lat**2) ** 1.5
    eta2 = nu / rho - 1

    VII = tan_lat / (2 * rho * nu)
    VIII = tan_lat / (24 * rho * nu**3) * (5 + 3 * tan_lat**2 + eta2 - 9 * tan_lat**2 * eta2)
    IX = tan_lat / (720 * rho * nu**5) * (61 + 90 * tan_lat**2 + 45 * tan_lat**4)
    X = 1 / (cos_lat * nu)
    XI = 1 / (cos_lat * 6 * nu**3) * (nu / rho + 2 * tan_lat**2)
    XII = 1 / (cos_lat * 120 * nu**5) * (5 + 28 * tan_lat**2 + 24 * tan_lat**4)
    XIIA = 1 / (cos_lat * 5040 * nu**7) * (61 + 662 * tan_lat**2 + 1320 * tan_lat**4 + 720 * tan_lat**6)

    dE = easting - E0
    lat = lat - VII * dE**2 + VIII * dE**4 - IX * dE**6
    lon = lon0 + X * dE - XI * dE**3 + XII * dE**5 - XIIA * dE**7

    # Helmert transform from OSGB36 to WGS84
    lat_deg = math.degrees(lat)
    lon_deg = math.degrees(lon)

    # Approximate correction (good enough for mapping)
    # Full Helmert would use 7-parameter transform, but this shortcut is ~5m accurate
    sin_lat2 = math.sin(lat)
    cos_lat2 = math.cos(lat)
    sin_lon = math.sin(lon)
    cos_lon = math.cos(lon)

    # Helmert parameters OSGB36 -> WGS84
    tx, ty, tz = 446.448, -125.157, 542.060
    s = -20.4894e-6
    rx = math.radians(0.1502 / 3600)
    ry = math.radians(0.2470 / 3600)
    rz = math.radians(0.8421 / 3600)

    nu_h = a / math.sqrt(1 - e2 * sin_lat2**2)
    x1 = (nu_h + 0) * cos_lat2 * cos_lon
    y1 = (nu_h + 0) * cos_lat2 * sin_lon
    z1 = (nu_h * (1 - e2) + 0) * sin_lat2

    x2 = tx + (1 + s) * x1 + (-rz) * y1 + (ry) * z1
    y2 = ty + (rz) * x1 + (1 + s) * y1 + (-rx) * z1
    z2 = tz + (-ry) * x1 + (rx) * y1 + (1 + s) * z1

    # WGS84 ellipsoid
    a2, b2 = 6378137.0, 6356752.3142
    e2_wgs = 1 - (b2 * b2) / (a2 * a2)
    p = math.sqrt(x2**2 + y2**2)
    wgs_lat = math.atan2(z2, p * (1 - e2_wgs))
    for _ in range(10):
        nu_wgs = a2 / math.sqrt(1 - e2_wgs * math.sin(wgs_lat) ** 2)
        wgs_lat_new = math.atan2(z2 + e2_wgs * nu_wgs * math.sin(wgs_lat), p)
        if abs(wgs_lat_new - wgs_lat) < 1e-12:
            break
        wgs_lat = wgs_lat_new

    wgs_lon = math.atan2(y2, x2)
    return round(math.degrees(wgs_lat), 6), round(math.degrees(wgs_lon), 6)


def extract_paths(filepath):
    with open(filepath, "r", encoding="utf-8", errors="replace") as f:
        content = f.read()

    features = []

    # Match each path block: path_id, path_name, path_distance, then mf_s() coordinate calls
    path_pattern = re.compile(
        r'\.path_id\s*=\s*(\d+);\s*'
        r".*?\.path_name\s*=\s*'(.*?)';\s*"
        r".*?\.path_name_url\s*=\s*'(.*?)';\s*"
        r".*?\.path_distance\s*=\s*'(.*?)';\s*"
        r"(.*?)var lineString",
        re.DOTALL,
    )

    coord_pattern = re.compile(r"mf_s\((\d+),(\d+)\)")

    for match in path_pattern.finditer(content):
        path_id = int(match.group(1))
        path_name = html.unescape(match.group(2))
        path_name_url = match.group(3)
        path_distance_raw = match.group(4).replace("&nbsp;", " ").strip()
        coord_block = match.group(5)

        coords = []
        for cm in coord_pattern.finditer(coord_block):
            easting, northing = int(cm.group(1)), int(cm.group(2))
            lat, lng = osgb36_to_wgs84(easting, northing)
            coords.append([lng, lat])  # GeoJSON is [lng, lat]

        if not coords:
            continue

        feature = {
            "type": "Feature",
            "properties": {
                "id": path_id,
                "name": path_name,
                "name_url": path_name_url,
                "distance": path_distance_raw,
            },
            "geometry": {
                "type": "LineString",
                "coordinates": coords,
            },
        }
        features.append(feature)

    return {
        "type": "FeatureCollection",
        "features": features,
    }


if __name__ == "__main__":
    geojson = extract_paths("search_by_path.php")
    output_file = "ldwa_paths.geojson"
    with open(output_file, "w") as f:
        json.dump(geojson, f)
    print(f"Extracted {len(geojson['features'])} paths to {output_file}")
