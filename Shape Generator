import ezdxf
import numpy as np
import scipy.spatial as ss
import random
import shapely.geometry as sg
import shapely.affinity as sa
import shapely.ops as so
import shapely.speedups
shapely.speedups.disable()
import scipy.stats as sst

def hexagonal_packing_coordinates(rows, cols, radius=1, separation=0):
    diam = radius * 2
    sep_factor = diam + separation

    x, y = np.meshgrid(np.arange(-cols // 2, cols // 2 + 1) * sep_factor, 
                       np.arange(-rows // 2, rows // 2 + 1) * sep_factor * np.sin(np.pi / 3))
    
    x[::2] += radius

    centers = list(zip(x.ravel(), y.ravel()))
    return centers

def generate_hexagon_vertices(area, centers):
    s = np.sqrt(2 * area / (3 * np.sqrt(3)))
    R = s / np.sqrt(3)

    vertices = []
    for i in range(6):
        angle = i * np.pi / 3
        x = R * np.sin(angle)
        y = R * np.cos(angle)
        vertices.append((x, y))
    
    hexagons = []
    for center in centers:
        hexagon = [(x + center[0], y + center[1]) for x, y in vertices]
        hexagons.append(hexagon)
    
    return hexagons

def scramble_centers(centers, area, scramble_factor=0.05, seed=69):
    random.seed(seed)
    max_deviation = area * scramble_factor
    
    scrambled_centers = []
    for coord in centers:
        x = coord[0] + random.uniform(-max_deviation, max_deviation)
        y = coord[1] + random.uniform(-max_deviation, max_deviation)
        scrambled_centers.append((x, y))
    
    return scrambled_centers

def generate_voronoi_vertices(centers):
    vor = ss.Voronoi(centers)
    vertices = vor.vertices
    regions = vor.regions
    finite_regions = [region for region in regions if -1 not in region]
    polygons = [vertices[region] for region in finite_regions]
    return polygons

def uniform_voronoi_area(polygons, target_area, separation):
    scaled_polygons = []

    for polygon in polygons:
        if polygon.shape[0] < 3:
            continue
        shapely_poly = sg.Polygon(polygon)
        scaled_poly = sa.scale(shapely_poly, np.sqrt(target_area / shapely_poly.area), np.sqrt(target_area / shapely_poly.area), origin="centroid")
        scaled_polygons.append(scaled_poly)
    
    while True:
        overlap = False
        for i, poly in enumerate(scaled_polygons):
            for j, other_poly in enumerate(scaled_polygons[i+1:]):
                if poly.intersects(other_poly):
                    overlap = True
                    distance = separation - poly.distance(other_poly)

                    vec = sg.Point(other_poly.centroid.x - poly.centroid.x, other_poly.centroid.y - poly.centroid.y)
                    tolerance = 1e-8  # a small tolerance value
                    vec_length = vec.length + tolerance
                    vec_x, vec_y = vec.x / vec_length * distance / 2, vec.y / vec_length * distance / 2
                    vec = sg.Point(vec_x, vec_y)

                    scaled_polygons[i] = sa.translate(scaled_polygons[i], vec.x, vec.y)
                    scaled_polygons[j+i+1] = sa.translate(scaled_polygons[j+i+1], -vec.x, -vec.y)
        if not overlap:
            break
    
    final_out = [np.array(list(p.exterior.coords)) for p in scaled_polygons]

    return final_out

def normal_dist_point_generator(width, height, std, num):
    # Define the mean and standard deviation of the normal distribution
    mean = [width/2, height/2]
    std = 10

    # Generate random points following a normal distribution
    x = np.random.normal(mean[0], std, num)
    y = np.random.normal(mean[1], std, num)

    # Create a 2D array of the points
    points = np.column_stack((x, y))
    return points

def generate_voronoi_vertices(centers):
    vor = ss.Voronoi(centers)
    vertices = vor.vertices
    regions = vor.regions
    finite_regions = [region for region in regions if -1 not in region]
    polygons = [vertices[region] for region in finite_regions]
    return polygons

def shrink_polygons(points, polygons, shrink):
    # for i, center in enumerate(points):
    #     for j, vertex in enumerate(vertices):
    #         dist = ss.distance.euclidean(center, vertex)
    #         new_vertex = vertex - ((vertex - center) / dist) * shrink
    #         vertices[j] = new_vertex
    scaled_polygons = []

    for polygon in polygons:
        # for point in points:
        if polygon.shape[0] < 3:
            continue

        shapely_poly = sg.Polygon(polygon)
        # shapely_point = sg.Point(point)
        scaled_poly = sa.scale(shapely_poly, shrink , shrink , origin="centroid")
        scaled_polygons.append(scaled_poly)
    
    final_out = [np.array(list(p.exterior.coords)) for p in scaled_polygons]

    return final_out
    # return scaled_polygons

def extr_dist_point_generator(width, height, std, num):
    # Define the mean and standard deviation of the normal distribution
    mean = [width/2, height/2]
    std = 10

    # Generate random points following a normal distribution
    x = np.random.normal(mean[0], std, num)
    y = x + np.random.normal(0, std/2, num)

    x = (x - np.min(x)) / (np.max(x) - np.min(x)) * height
    y = (y - np.min(y)) / (np.max(y) - np.min(y)) * width

    # Create a 2D array of the points
    points = np.column_stack((x, y))
    return points

def generate_voronoi_vertices(centers):
    vor = ss.Voronoi(centers)
    vertices = vor.vertices
    regions = vor.regions
    finite_regions = [region for region in regions if -1 not in region]
    polygons = [vertices[region] for region in finite_regions]
    return polygons

def shrink_polygons(points, polygons, shrink):
    scaled_polygons = []

    for polygon in polygons:
        # for point in points:
        if polygon.shape[0] < 3:
            continue

        shapely_poly = sg.Polygon(polygon)
        scaled_poly = sa.scale(shapely_poly, shrink , shrink , origin="centroid")
        scaled_polygons.append(scaled_poly)
    
    final_out = [np.array(list(p.exterior.coords)) for p in scaled_polygons]

    return final_out

# The list of adjustable parameters were as follows:

# Parameters
area = 50
separation = 0.2
no_of_rows = 50
no_of_cols = 50
random_seed = 100 # do not change this unless you want a different generation
scramble_factor = 0.2

width = 45 * 2 * 3
height = 55 * 2 * 3
stdev = 10
num_of_points = 80
shrink = 0.8

# To generate the DXF files with shape arrays:

# =======================
# Circle Array
# =======================

# Create a document
circle_doc = ezdxf.new('R2018') 
circle_msp = circle_doc.modelspace()

# Derive radius from cell area
circle_radius = np.sqrt(area / np.pi)

# Generate circle center coordinates
circle_centers = hexagonal_packing_coordinates(rows=no_of_rows, cols=no_of_cols, radius = circle_radius, separation=separation) # inclusion of separation distance for centerpoints

# Add to DXF file
for coord in circle_centers:
    circle_msp.add_circle(center=coord, radius=circle_radius)

# Saves document
circle_doc.saveas("circle.dxf")

# =======================
# Hexagonal/Honeycomb Array
# =======================

# Create a document
hexagon_doc = ezdxf.new('R2018') 
hexagon_msp = hexagon_doc.modelspace()

# Calculate side lengths
s = np.sqrt(2 * area / (3 * np.sqrt(3)))
hex_R = s / np.sqrt(3)

# Generate hexagon coordinates 
hex_array_coords = hexagonal_packing_coordinates(rows=no_of_rows, cols=no_of_cols, radius = hex_R, separation=separation) # inclusion of separation distance for centerpoints
hexagons = generate_hexagon_vertices(area, hex_array_coords)

# Add to DXF file
for hex in hexagons:
    hexagon_msp.add_polyline2d(hex, close = True)

# Adds weight point
for point in hex_array_coords:
    hexagon_msp.add_point(point)

# Saves document
hexagon_doc.saveas("honeycomb.dxf")

# =======================
# Voronoi Array
# =======================

# Create a document
voro_doc = ezdxf.new('R2018') 
voro_msp = voro_doc.modelspace()

# Derive radius from cell area
circle_radius = np.sqrt(area / np.pi) # to minimise the "ultra random lul" behaviour, start with a regular array

# Generate hexagonal packing  
voro_array_coords = hexagonal_packing_coordinates(rows=no_of_rows, cols=no_of_cols, radius = circle_radius, separation=separation) # inclusion of separation distance for centerpoints
voro_scrambled = scramble_centers(centers=voro_array_coords, area=area, scramble_factor=scramble_factor, seed=random_seed)
voro_polygons = generate_voronoi_vertices(voro_scrambled)
voro_uniform_area_polys = uniform_voronoi_area(voro_polygons, area, separation)

# Add to DXF file
for poly in voro_uniform_area_polys:
    voro_msp.add_polyline2d(poly, close = True)

# Adds weight point
for point in voro_scrambled:
    voro_msp.add_point(point)

# Saves document
voro_doc.saveas("voronois.dxf")

# =======================
# Radial Normal Dist Array
# =======================

# Create a document
norm_doc = ezdxf.new('R2018') 
norm_msp = norm_doc.modelspace()

# Derive radius from cell area
circle_radius = np.sqrt(area / np.pi) # to minimise the "ultra random lul" behaviour, start with a regular array

# Generate hexagonal packing  

norm_coords = normal_dist_point_generator(width, height, stdev, num_of_points)
norm_polygons = generate_voronoi_vertices(norm_coords)
norm_polygons_shrunk = shrink_polygons(norm_coords, norm_polygons, shrink)
# norm_polygons_sep = uniform_voronoi_area(norm_polygons, area, separation)

# Add to DXF file
for poly in norm_polygons_shrunk:
    norm_msp.add_polyline2d(poly, close = True)

# Adds weight point
for point in norm_coords:
    norm_msp.add_point(point)

# Saves document
norm_doc.saveas("normonois new seed .dxf")

# =======================
# Extrusion Normal Dist Array
# =======================

# Create a document
extr_doc = ezdxf.new('R2018') 
extr_msp = extr_doc.modelspace()

# Derive radius from cell area
circle_radius = np.sqrt(area / np.pi) # to minimise the "ultra random lul" behaviour, start with a regular array

# Generate hexagonal packing  

extr_coords = extr_dist_point_generator(width, height, stdev, num_of_points)
extr_polygons = generate_voronoi_vertices(extr_coords)
extr_polygons_shrunk = shrink_polygons(extr_coords, extr_polygons, shrink)

# Add to DXF file
for poly in extr_polygons_shrunk:
    extr_msp.add_polyline2d(poly, close = True)

# Adds weight point
for point in extr_coords:
    extr_msp.add_point(point)

# Saves document
extr_doc.saveas("extruderois.dxf")
