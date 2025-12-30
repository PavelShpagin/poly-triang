"""
Fix for the compute_left_edges bug in triangulation.py

The issue: When finding left_edge for a vertex, incoming edges (edges that END
at the vertex) should be removed BEFORE querying, not after.

This script patches the compute_left_edges method.
"""

import os

# Read the current file
with open("note/triangulation.py", "r") as f:
    content = f.read()

# The buggy section (finds left edge BEFORE removing incoming edge):
old_code = '''        # Find left edge BEFORE modifying active edges
        left_e = find_left_edge(v.x, y)
        if left_e is not None:
            self.left_edge[idx] = left_e
        
        # Update active edges
        if incoming_edge is not None:
            remove_edge(incoming_edge)
        if outgoing_edge is not None:
            insert_edge(outgoing_edge, y)'''

# Fixed version (removes incoming edge BEFORE finding left edge):
new_code = '''        # IMPORTANT: Remove incoming edge BEFORE finding left edge
        # Otherwise, edges ending at this vertex may be incorrectly selected
        if incoming_edge is not None:
            remove_edge(incoming_edge)
        
        # Now find the left edge
        left_e = find_left_edge(v.x, y)
        if left_e is not None:
            self.left_edge[idx] = left_e
        
        # Insert outgoing edge
        if outgoing_edge is not None:
            insert_edge(outgoing_edge, y)'''

if old_code in content:
    content = content.replace(old_code, new_code)
    with open("note/triangulation.py", "w") as f:
        f.write(content)
    print("Patched successfully!")
else:
    print("Could not find the code to patch. Manual fix needed.")
    print("\nLooking for similar patterns...")
    # Try to find where the issue might be
    if "Find left edge BEFORE" in content:
        print("Found 'Find left edge BEFORE' comment")
    if "find_left_edge" in content:
        print("Found 'find_left_edge' function call")

