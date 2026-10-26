import pya
import math


def create_triangle(cell, layer, cx, cy, a1, a2, theta1):
    """Generates a triangular shape in a KLayout layout with the given center, side lengths, and angles."""
    
    # Convert angles from degrees to radians
    theta1_rad = math.radians(theta1)
    
    
    
    p1 = pya.DPoint(cx, cy)
    p2 = pya.DPoint(cx - a1, cy)
    p3 = pya.DPoint(p1.x - a2 * math.cos(theta1_rad) , cy - a2 * math.sin(theta1_rad))

    # Define the polygon as the triangle
    polygon = pya.DPolygon([p1, p2, p3])
    cell.shapes(layer).insert(polygon)

    print(f"Triangle generated at center: ({cx}, {cy})")


