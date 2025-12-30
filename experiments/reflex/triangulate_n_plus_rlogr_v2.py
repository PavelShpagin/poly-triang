import math
import random
from dataclasses import dataclass
from enum import Enum
from typing import Dict, List, Tuple


class VertexType(Enum):
    START = 1
    END = 2
    SPLIT = 3
    MERGE = 4
    REGULAR_LEFT = 5
    REGULAR_RIGHT = 6


@dataclass
class Vertex:
    x: float
    y: float
    idx: int
    type: VertexType = None
    prev_idx: int = -1
    next_idx: int = -1


@dataclass
class Edge:
    upper_idx: int
    lower_idx: int
    x_at_y: float
    helper_idx: int = -1
    prev: "Edge" = None
    next: "Edge" = None
    in_tree: bool = False


def cross(o, a, b):
    return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x)


def x_intercept(v1, v2, y):
    if abs(v1.y - v2.y) < 1e-12:
        return (v1.x + v2.x) / 2
    t = (y - v1.y) / (v2.y - v1.y)
    return v1.x + t * (v2.x - v1.x)


class SkipNode:
    __slots__ = ("key", "edge", "forward")

    def __init__(self, key, edge, level):
        self.key = key
        self.edge = edge
        self.forward = [None] * (level + 1)


class SkipList:
    MAX_LEVEL = 16

    def __init__(self):
        self.header = SkipNode(float("-inf"), None, self.MAX_LEVEL)
        self.level = 0
        self.size = 0

    def _random_level(self):
        lvl = 0
        while random.random() < 0.5 and lvl < self.MAX_LEVEL:
            lvl += 1
        return lvl

    def insert(self, key, edge):
        update = [None] * (self.MAX_LEVEL + 1)
        cur = self.header
        for i in range(self.level, -1, -1):
            while cur.forward[i] and cur.forward[i].key < key:
                cur = cur.forward[i]
            update[i] = cur
        lvl = self._random_level()
        if lvl > self.level:
            for i in range(self.level + 1, lvl + 1):
                update[i] = self.header
            self.level = lvl
        node = SkipNode(key, edge, lvl)
        for i in range(lvl + 1):
            node.forward[i] = update[i].forward[i]
            update[i].forward[i] = node
        self.size += 1
        edge.in_tree = True

    def remove(self, key):
        update = [None] * (self.MAX_LEVEL + 1)
        cur = self.header
        for i in range(self.level, -1, -1):
            while cur.forward[i] and cur.forward[i].key < key:
                cur = cur.forward[i]
            update[i] = cur
        target = cur.forward[0]
        if target and abs(target.key - key) < 1e-9:
            for i in range(self.level + 1):
                if update[i].forward[i] != target:
                    break
                update[i].forward[i] = target.forward[i]
            self.size -= 1

    def predecessor(self, key):
        cur = self.header
        for i in range(self.level, -1, -1):
            while cur.forward[i] and cur.forward[i].key < key:
                cur = cur.forward[i]
        return cur.edge if cur != self.header else None


class EdgeList:
    def __init__(self):
        self.head = Edge(-1, -1, float("-inf"))
        self.tail = Edge(-1, -1, float("inf"))
        self.head.next = self.tail
        self.tail.prev = self.head

    def insert_after(self, pos, edge):
        edge.prev = pos
        edge.next = pos.next
        pos.next.prev = edge
        pos.next = edge

    def remove(self, edge):
        edge.prev.next = edge.next
        edge.next.prev = edge.prev


class Triangulator:
    def __init__(self, vertices: List[Vertex]):
        self.V = vertices
        self.n = len(vertices)
        self.L = EdgeList()
        self.T = SkipList()
        self.diagonals = []
        self.r = 0
        self.edge_of = {}
        self.left_edge_of = {}

    def classify_vertices(self):
        for i in range(self.n):
            v = self.V[i]
            v.prev_idx = (i - 1) % self.n
            v.next_idx = (i + 1) % self.n
            prev, nxt = self.V[v.prev_idx], self.V[v.next_idx]
            below_prev = prev.y < v.y
            below_next = nxt.y < v.y
            reflex = cross(prev, v, nxt) < 0

            if below_prev and below_next:
                v.type = VertexType.SPLIT if reflex else VertexType.START
                if reflex:
                    self.r += 1
            elif not below_prev and not below_next:
                v.type = VertexType.MERGE if reflex else VertexType.END
                if reflex:
                    self.r += 1
            else:
                v.type = (
                    VertexType.REGULAR_LEFT
                    if prev.y > nxt.y
                    else VertexType.REGULAR_RIGHT
                )

    def add_landmark(self, edge):
        if edge and not edge.in_tree and edge != self.L.head and edge != self.L.tail:
            self.T.insert(edge.x_at_y, edge)

    def find_left_edge(self, x):
        landmark = self.T.predecessor(x)
        cur = landmark if landmark else self.L.head
        while cur.next != self.L.tail and cur.next.x_at_y < x:
            if cur.next != self.L.tail:
                self.add_landmark(cur.next)
            cur = cur.next
        return cur if cur != self.L.head else None

    def handle_start(self, v, y):
        nxt = self.V[v.next_idx]
        x = x_intercept(v, nxt, y)
        edge = Edge(v.idx, v.next_idx, x, helper_idx=v.idx)
        left = self.find_left_edge(x)
        pos = left if left else self.L.head
        self.L.insert_after(pos, edge)
        self.add_landmark(edge)
        self.edge_of[v.next_idx] = edge
        self.left_edge_of[v.next_idx] = edge

    def handle_end(self, v, y):
        edge = self.edge_of.get(v.idx)
        if edge:
            if edge.in_tree:
                self.T.remove(edge.x_at_y)
            self.L.remove(edge)
            del self.edge_of[v.idx]

    def handle_split(self, v, y):
        left = self.find_left_edge(v.x)
        if left and left != self.L.head:
            self.add_landmark(left)
            if left.helper_idx >= 0:
                self.diagonals.append((v.idx, left.helper_idx))
            left.helper_idx = v.idx

        nxt = self.V[v.next_idx]
        x = x_intercept(v, nxt, y)
        edge = Edge(v.idx, v.next_idx, x, helper_idx=v.idx)
        pos = left if left else self.L.head
        self.L.insert_after(pos, edge)
        self.edge_of[v.next_idx] = edge
        self.left_edge_of[v.next_idx] = edge

    def handle_merge(self, v, y):
        edge_in = self.edge_of.get(v.idx)
        if edge_in:
            if edge_in.helper_idx >= 0:
                helper = self.V[edge_in.helper_idx]
                if helper.type == VertexType.MERGE:
                    self.diagonals.append((v.idx, edge_in.helper_idx))
            if edge_in.in_tree:
                self.T.remove(edge_in.x_at_y)
            self.L.remove(edge_in)
            del self.edge_of[v.idx]

        left = self.find_left_edge(v.x)
        if left and left != self.L.head:
            self.add_landmark(left)
            if left.helper_idx >= 0:
                helper = self.V[left.helper_idx]
                if helper.type == VertexType.MERGE:
                    self.diagonals.append((v.idx, left.helper_idx))
            left.helper_idx = v.idx

    def handle_regular_left(self, v, y):
        edge_in = self.edge_of.get(v.idx)
        if edge_in:
            if edge_in.helper_idx >= 0:
                helper = self.V[edge_in.helper_idx]
                if helper.type == VertexType.MERGE:
                    self.diagonals.append((v.idx, edge_in.helper_idx))

            was_landmark = edge_in.in_tree
            prev_pos = edge_in.prev
            if was_landmark:
                self.T.remove(edge_in.x_at_y)
            self.L.remove(edge_in)
            del self.edge_of[v.idx]

            nxt = self.V[v.next_idx]
            x = x_intercept(v, nxt, y)
            edge = Edge(v.idx, v.next_idx, x, helper_idx=v.idx)
            self.L.insert_after(prev_pos, edge)
            self.edge_of[v.next_idx] = edge
            if was_landmark:
                self.T.insert(x, edge)

    def handle_regular_right(self, v, y):
        left = self.left_edge_of.get(v.idx)
        if left and left != self.L.head:
            if left.helper_idx >= 0:
                helper = self.V[left.helper_idx]
                if helper.type == VertexType.MERGE:
                    self.diagonals.append((v.idx, left.helper_idx))
            left.helper_idx = v.idx
            nxt = self.V[v.next_idx]
            if nxt.type == VertexType.REGULAR_RIGHT:
                self.left_edge_of[v.next_idx] = left

    def triangulate(self) -> List[Tuple[int, int]]:
        self.classify_vertices()
        order = sorted(range(self.n), key=lambda i: -self.V[i].y)

        for idx in order:
            v = self.V[idx]
            y = v.y - 1e-9
            handlers = {
                VertexType.START: self.handle_start,
                VertexType.END: self.handle_end,
                VertexType.SPLIT: self.handle_split,
                VertexType.MERGE: self.handle_merge,
                VertexType.REGULAR_LEFT: self.handle_regular_left,
                VertexType.REGULAR_RIGHT: self.handle_regular_right,
            }
            handlers[v.type](v, y)

        return self.diagonals


def star_polygon(n):
    vertices = []
    for i in range(n):
        angle = 2 * math.pi * i / n - math.pi / 2
        radius = 2 if i % 2 == 0 else 1
        vertices.append(
            Vertex(radius * math.cos(angle), radius * math.sin(angle), i)
        )
    return vertices


def smooth_polygon(n):
    vertices = []
    for i in range(n):
        angle = 2 * math.pi * i / n
        radius = 1 + 0.05 * math.sin(3 * angle)
        vertices.append(
            Vertex(radius * math.cos(angle), radius * math.sin(angle), i)
        )
    return vertices


def benchmark():
    print("Polygon Triangulation: O(n + r log r)")
    print("=" * 60)
    print(f"{'n':>8} {'r':>6} {'diagonals':>10}")
    print("-" * 60)

    for n in [100, 500, 1000, 2000, 5000, 10000]:
        verts = star_polygon(n)
        tri = Triangulator(verts)
        diags = tri.triangulate()
        print(f"{n:>8} {tri.r:>6} {len(diags):>10}")

    print()
    for n in [100, 500, 1000, 2000, 5000, 10000]:
        verts = smooth_polygon(n)
        tri = Triangulator(verts)
        diags = tri.triangulate()
        print(f"{n:>8} {tri.r:>6} {len(diags):>10}")


if __name__ == "__main__":
    benchmark()
