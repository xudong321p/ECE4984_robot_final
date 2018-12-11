#!/usr/bin/env python

import copy
import math

"""
the graph object
"""


class Graph:
    def __init__(self, _all_vertex, _start, _end):
        self.start = _start
        self.end = _end
        self.all_vertex = _all_vertex
        self.edges = [[-1 for i in range(self.all_vertex)] for i in range(self.all_vertex)]
        self.closed_vertex = dict()

    def reset_state(self):
        self.closed_vertex = dict()

    def set_coords(self, _coords):
        self.coords = copy.deepcopy(_coords)

    def put_into_closed_vertex(self, n, distance):
        self.closed_vertex[n] = distance

    """
    get point connected to the n
    """

    def get_distance_direct_to(self, n):
        col = n
        dist = dict()
        for i in range(self.all_vertex):
            if self.edges[i][col] != -1 and not self.closed_vertex.has_key(i):
                dist[i] = self.edges[i][col]
        return dist

    """
    get point the n connect to
    """

    def get_distance_direct_from(self, n):
        row = n
        dist = dict()
        for i in range(self.all_vertex):
            if self.edges[row][i] != -1 and not self.closed_vertex.has_key(i):
                dist[i] = self.edges[row][i]
        return dist

    """
    calculate the Eucliden distance
    """

    def get_Euclidean_distance(self, node1, node2):
        coord1 = self.coords[node1]
        coord2 = self.coords[node2]
        dx = coord1[0] - coord2[0]
        dy = coord1[1] - coord2[1]
        d = math.sqrt(dx * dx + dy * dy)
        return d


"""
the temporary point record
"""


class Record:
    def __init__(self, _from, _to, _dist):
        self._from = _from
        self.to = _to
        self.dist = _dist
        self.F = -1


"""
read from the file
"""


def handle_input(index):
    file_name = "input_" + str(index) + ".txt"
    coord_name = "coords_" + str(index) + ".txt"

    graph = read_from_file(file_name)
    coords = read_from_coords(coord_name)
    graph.set_coords(coords)
    return graph


def read_from_file(file_name):
    fp = open(file_name, "r")
    line = fp.readline().strip()
    all_vertex = int(line, 10)
    line = fp.readline().strip()
    start = int(line, 10)
    line = fp.readline().strip()
    end = int(line, 10)

    graph = Graph(all_vertex, start - 1, end - 1)

    for line in fp:
        line = line.strip()
        tmp_array = line.split()
        if len(tmp_array) < 3:
            continue
        _from = int(tmp_array[0], 10) - 1
        _to = int(tmp_array[1], 10) - 1
        _w = float(tmp_array[2])
        graph.edges[_from][_to] = _w

    fp.close()
    # print_graph(graph)
    return graph


def read_from_coords(file_name):
    fp = open(file_name)
    coords = dict()
    index = 0
    for line in fp:
        line = line.strip()
        tmp = line.split()
        pos = [float(tmp[0]), float(tmp[1])]
        coords[index] = pos
        index += 1
    return coords


def print_graph(_graph):
    print("----number of vertex----%d----" % (_graph.all_vertex))
    print("----start----%d----" % (_graph.start))
    print("----end----%d----" % (_graph.end))
    print("----weight--------------------")
    for i in range(_graph.all_vertex):
        for j in range(_graph.all_vertex):
            if _graph.edges[i][j] != -1:
                print
                i, j, _graph.edges[i][j]


"""
dijkstra method
"""


def dijkstra(_graph, _start, _end):
    if _start == _end:
        return [[_start, _end], 0]

    current = _start
    open_vertex = dict()
    record_list = dict()
    for i in range(_graph.all_vertex):
        record_list[i] = Record(i, i, -1)

    _graph.closed_vertex[_start] = 0
    record_list[_start]._from = _start
    record_list[_start].to = _start
    record_list[_start].dist = 0
    while current != _end:
        vertex_list = _graph.get_distance_direct_from(current)
        current_dist = record_list[current].dist
        for vertex in vertex_list.keys():
            if open_vertex.has_key(vertex):
                if record_list[vertex].dist > current_dist + vertex_list[vertex]:
                    open_vertex[vertex] = current_dist + vertex_list[vertex]
                    record_list[vertex].dist = current_dist + vertex_list[vertex]
                    record_list[vertex]._from = current
            else:
                open_vertex[vertex] = current_dist + vertex_list[vertex]
                record_list[vertex].dist = current_dist + vertex_list[vertex]
                record_list[vertex]._from = current
        # new current
        tmp_dist = -1
        for vertex in open_vertex.keys():
            if tmp_dist == -1:
                tmp_dist = open_vertex[vertex]
                tmp_vertex = vertex
            elif tmp_dist > open_vertex[vertex]:
                tmp_dist = open_vertex[vertex]
                tmp_vertex = vertex
        if tmp_dist != -1:
            current = tmp_vertex
            _graph.put_into_closed_vertex(current, tmp_dist)
            del open_vertex[current]
        else:
            break

    route_list = get_route(record_list, _start, _end)
    dist = _graph.closed_vertex[_end]
    return [route_list, dist]


"""
A* method
"""


def AStar(_graph, _start, _end):
    current = _start
    open_vertex = dict()
    record_list = dict()
    for i in range(_graph.all_vertex):
        record_list[i] = Record(i, i, -1)

    _graph.closed_vertex[_start] = 0
    record_list[_start]._from = _start
    record_list[_start].to = _start
    record_list[_start].dist = 0
    record_list[_start].F = 0
    while current != _end:
        vertex_list = _graph.get_distance_direct_from(current)
        current_dist = record_list[current].dist
        for vertex in vertex_list.keys():
            G = current_dist + vertex_list[vertex]
            F = G + _graph.get_Euclidean_distance(vertex, _end)
            if open_vertex.has_key(vertex):
                if open_vertex[vertex] > F:
                    open_vertex[vertex] = F
                    record_list[vertex].dist = G
                    record_list[vertex].F = F
                    record_list[vertex]._from = current
            else:
                open_vertex[vertex] = F
                record_list[vertex].dist = G
                record_list[vertex].F = F
                record_list[vertex]._from = current
        # new current
        tmp_dist = -1
        for vertex in open_vertex.keys():
            if tmp_dist == -1:
                tmp_dist = open_vertex[vertex]
                tmp_vertex = vertex
            elif tmp_dist > open_vertex[vertex]:
                tmp_dist = open_vertex[vertex]
                tmp_vertex = vertex
        if tmp_dist != -1:
            current = tmp_vertex
            tmp_dist = record_list[current].dist
            _graph.put_into_closed_vertex(current, tmp_dist)
            del open_vertex[current]
        else:
            break

    route_list = get_route(record_list, _start, _end)
    dist = _graph.closed_vertex[_end]
    return [route_list, dist]


"""
floyd method
"""


def floyd(_graph, _start, _end):
    distance = [[-1 for i in range(_graph.all_vertex)] for i in range(_graph.all_vertex)]
    path = [[-1 for i in range(_graph.all_vertex)] for i in range(_graph.all_vertex)]
    for row in range(_graph.all_vertex):
        for col in range(_graph.all_vertex):
            if row == col:
                distance[row][col] = 0
                path[row][col] = row
            else:
                distance[row][col] = _graph.edges[row][col]
                if _graph.edges[row][col] == -1:
                    path[row][col] = -1
                else:
                    path[row][col] = row

    for k in range(_graph.all_vertex):
        for row in range(_graph.all_vertex):
            for col in range(_graph.all_vertex):
                if distance[row][k] != -1 and distance[k][col] != -1:
                    dist = distance[row][k] + distance[k][col]
                    if distance[row][col] == -1 or distance[row][col] > dist:
                        distance[row][col] = dist
                        path[row][col] = k

    route_list = []
    end = _end
    while end != _start:
        route_list.append(end)
        end = path[_start][end]
    route_list.append(_start)
    route_list.reverse()

    dist = distance[_start][_end]
    return [route_list, dist]


"""
find the path
"""


def get_route(record_list, _start, _end):
    current = _end
    route_list = []
    while current != _start:
        route_list.append(current)
        current = record_list[current]._from
    route_list.append(_start)
    route_list.reverse()
    return route_list


"""
write to the file
"""


def handle_output(dist, iter_nums, counts):
    fp = open("output_costs.txt", "w")
    for i in range(counts):
        fp.write(str(dist[i][0]) + " " + str(dist[i][1]) + " " + str(dist[i][2]))
        fp.write("\n")
    fp.close()

    fp = open("output_numiters.txt", "w")
    for i in range(counts):
        fp.write(str(iter_nums[i][0]) + " " + str(iter_nums[i][1]) + " " + str(iter_nums[i][2]))
        fp.write("\n")
    fp.close()


def main():
    counts = 3
    dist = [[0 for i in range(3)] for i in range(counts)]
    iter_nums = [[0 for i in range(3)] for i in range(counts)]

    for index in range(counts):
        graph = handle_input(index + 1)
        [route_list, _dist] = dijkstra(graph, graph.start, graph.end)
        dist[index][0] = _dist
        iter_nums[index][0] = len(graph.closed_vertex.keys())
        del graph

        graph = handle_input(index + 1)
        [route_list, _dist] = AStar(graph, graph.start, graph.end)
        dist[index][1] = _dist
        iter_nums[index][1] = len(graph.closed_vertex.keys())
        del graph

        graph = handle_input(index + 1)
        [route_list, _dist] = floyd(graph, graph.start, graph.end)
        dist[index][2] = _dist
        iter_nums[index][2] = pow(graph.all_vertex, 3)
        del graph

    handle_output(dist, iter_nums, counts)


if __name__ == "__main__":
    main()
