'''
This code enables the visualization of (hopefully) non-trivial cycles of non-k-equal configurations of stars 
(see https://en.wikipedia.org/wiki/Star_(graph_theory)) using Tkinter.
                                                                                   . 
For example, this is a 3-star (also known a Y-graph, because it looks like a Y):   |  
                                                                                   .__.
                                                                                  /   
                                                                                 .
Here, we say that d=2.

For a given d-star, where d is the number of "arms" extending from the central vertex, we construct a visualization of the
aforementioned cycles by creating the Abrams-discretized model of the non-k-equal configuration space of n points (or "robots")
living on the graph.

For several reasons, this only really works 100% for 3-stars at the moment. 
TODO: [d>3] There is no single code to generate Abrams-discretized models for any given d-star, e.g. the code in conf_n_k_Y.py 
      only works for d=3. It would be great if the code in conf_n_k_Y.py could be generalized to build any d-star. That code would 
      easily be compatible with the code currently written in this file.
TODO: [d=2] The StarGraph class defined below doesn't work for d=2 (i.e. the graph that looks like a line/interval, or an I), because 
      that complex as defined in conf_n_k_I.py uses a vertex at one of the ENDS as its base point and creates "downstream" cubes from 
      that endpoint, whereas for d>2, "downstream" cubes are constructed based on the CENTER vertex. Videos for D_{3,3}I and D_{4,3}I
      were generated fairly easily by making some brute force changes to this code, but it would also be good for the StarGraph class
      to just support d=2 as is. 
'''
import math
import random
import time
import Tkinter as tk
from homology.conf_n_k_Y import the_complex as Y_COMPLEX
# from homology.conf_n_k_I import the_complex as I_COMPLEX

# We need sine and cosine functions to calculate the movement angle of robots moving along the graph.
sin = lambda degs: math.sin(math.radians(degs))
cos = lambda degs: math.cos(math.radians(degs))

class StarGraph(object):
    def __init__(self, starNum, n, k=2):
        '''
        :param starNum: The number of "arms" a given star has (see above description). If starNum = d then we call our graph a d-star.
        :type starNum: int
        :param n: An integer representing the number of moving points (or "robots") living in our graph Y. Must have n >= k.
        :type n: int
        :param k: The integer k for the non-k-equal configuration space. The default is 2 (no points can cross).
        :type k: int
        '''

        # Make sure the numbers are not wonky.
        if starNum < 2: raise ValueError(str(starNum) + "-stars don't really make sense. Make sure starNum is at least 2.")
        if starNum != 3: raise NotImplementedError("Sorry! The biggest stars implemented right now are 3-stars.")
        if k > n:       raise ValueError("The value k can only be at most as large as the value n.")

        self._FRAME_SIZE = 600.0 # Change this to suit the screen you are using.
        self._CENTER_POINT = self._FRAME_SIZE/2
        self._OFFSET_RADIUS = self._FRAME_SIZE/2 - 10

        # These keep track of which cycles is currently being displayed, and which position in the cycle the robots are at.
        self.current_cycle = 0
        self.current_position = 0

        self._edges = [[] for i in range(starNum)]
        cubical_complex = Y_COMPLEX(n,k) if starNum == 3 else I_COMPLEX(n,k)

        # TODO: Currently, when calling cubical_complex.sorted_n_cycles(n), n has to be manually changed to get the desired
        #       n-cycles. Perhaps there a should be a neater way to choose n, maybe based on the largest homology.
        self._cycles = cubical_complex.sorted_n_cycles(1)
        assert(len(self._cycles) > 0)

        # Helper setup methods
        self._tk_setup(starNum, n, k)
        self._init_points(starNum, n, k)
        self._init_robots(starNum, n, k)

        # Animate with Tkinter
        self._next_position(starNum, n, k)
        self._root.mainloop()

    def _tk_setup(self, starNum, n, k):
        '''
        Set up Tkinter.
        '''
        self._root = tk.Tk()
        self._root.title(str(starNum) + "-Star Graph with " + str(n) + " Robots " + "(Non-" + str(k) + "-Equal)")
        self.canvas = tk.Canvas(self._root, bg="white", height=self._FRAME_SIZE, width=self._FRAME_SIZE)
        next_cycle_button = tk.Button(self._root, text="Next", command=lambda: self._next_cycle_callback(starNum, n, k))
        next_cycle_button.pack()
        self.cycle_label = tk.Label(self._root, text="Cycle #"+str(self.current_cycle+1) + " of #"+str(len(self._cycles)) 
                                    + "\nPosition " + str(self.current_position+1) + " of " + str(len(self._cycles[self.current_cycle])))
        self.cycle_label.pack()

    def _init_points(self, starNum, n, k):
        '''
        Set up the the center vertex build edges from it to outer vertices.
        '''
        self.center = Vertex(self, self._CENTER_POINT, self._CENTER_POINT)
        num_edges = n - 1
        for i in range(num_edges):
            for j in range(starNum):
                angle = 360.0/starNum * (j+1)
                if i == 0:
                    start = self.center
                else:
                    current = self._edges[j]
                    start = current[len(current)-1].end
                end = Vertex(self, start.x+((self._OFFSET_RADIUS/num_edges)*sin(angle)), (start.y+(self._OFFSET_RADIUS/num_edges)*cos(angle)))
                edge = Edge(self, start, end)
                self._edges[j].append(edge)

    def _init_robots(self, starNum, n, k):
        '''
        Build the robots and put them at the start position.
        '''
        self.robots = []
        for i in range(n):
            point = self._get_point_for_robot(starNum, n, k, i)
            robot = Robot(self, point.x, point.y, i)
            if point != self.center:
                robot.centered = False
            self.robots.append(robot)
        self.canvas.pack()

    def _get_point_for_robot(self, starNum, n, k, i):
        '''
        At any given moment, get the point where Robot i is meant to be.
        '''
        dimension = len(self._cycles[0][0][1].tuple()) # TODO: is there a more elegant way to do this?
        num_robots = n
        positions = self._cycles[self.current_cycle][self.current_position][1].tuple()
        robot_position = positions[i*starNum:i*starNum+(dimension/n)]
        for j in range(len(robot_position)):
            if not robot_position[j] == (0, 0):
                startPosition = robot_position[j][0]
                endPosition   = robot_position[j][1]
                if startPosition == endPosition:
                    point = self._edges[j%starNum][endPosition-1].end
                else:
                    point = self._edges[j%starNum][endPosition-1].midpoint()
                return point
        return self.center

    def _move_robots(self, starNum, n, k):
        '''
        Helper method called by _next_position to move the robots.
        '''
        for i in range(len(self.robots)):
            point = self._get_point_for_robot(starNum, n, k, i)
            if self.robots[i].current_point.angle_to_point(self.center) == point.angle_to_point(self.center) or self.robots[i].current_point == self.center:
                if point == self.center:
                    self.robots[i].centered = True
                self.robots[i].move_to(point)
            else:
                self.robots[i].centered = False
                self.robots[i].move_to(self.center)

    def _next_position(self, starNum, n, k):
        '''
        This is the animation method called by Tkinter. It facilitates moving from one position in a cycle to the next.
        When it gets to the last position it simply moves back to the first position.
        '''
        self.current_position += 1
        position_changed = True
        for robot in self.robots:
            if robot.next_point is not None or (robot.current_point==self.center and not robot.centered):
                self.current_position -= 1
                position_changed = False
                break
        if self.current_position == len(self._cycles[self.current_cycle]):
            self.current_position = 0
        if position_changed:
            time.sleep(0.3)
            self.cycle_label['text'] = "Cycle #"+str(self.current_cycle+1) + " of #"+str(len(self._cycles)) + "\nPosition " + str(self.current_position+1) + " of " + str(len(self._cycles[self.current_cycle]))
            time.sleep(0.3)
        self._move_robots(starNum, n, k)
        self.animation_id = self.canvas.after(30, lambda: self._next_position(starNum, n, k))

    def _reset(self, starNum, n, k):
        '''
        This method is called every time the cycle is changed. It reconfigures the robots and runs the animation for the new cycle.
        '''
        self.canvas.after_cancel(self.animation_id)
        self.current_cycle += 1
        if self.current_cycle >= len(self._cycles)-1:
            self.current_cycle %= len(self._cycles)
        self.current_position = 0
        for robot in self.robots:
            self.canvas.delete(robot.oval)
        self._init_robots(starNum, n, k)
        self.current_position = -1 # A temporary fix
        self._next_position(starNum, n, k)

    def _next_cycle_callback(self, starNum, n, k):
        '''
        A callback method called every time the user clicks the "Next button on the Tkinter window.
        '''
        self._reset(starNum, n, k)
        self.cycle_label['text'] = "Cycle #"+str(self.current_cycle+1) + " of #"+str(len(self._cycles)) + "\nPosition " + str(self.current_position+1) + " of " + str(len(self._cycles[self.current_cycle]))

class Point(object):
    '''
    Point is a superclass for vertices and robots that defines their shared traits.
    '''
    # Frequently used sine and cosine values.
    _COS_0   = cos(0)
    _COS_180 = cos(180)
    _SIN_90  = sin(90)
    _SIN_270 = sin(270)

    # How large vertices and robots are.
    _POINT_RADIUS = 5

    def __init__(self, graph, x, y):
        self._graph = graph
        self._color = "black"

        # These values are important for robot movement:
        # - self.x and self.y encode the REAL coordinates of a point in the Tkinter window,
        # - self.current_point is a copy of this point that encodes the last position where this point was stationary.
        # - self.next_point is None if the robot is not currently going anywhere, and contains a point if it has a destination.
        self.x, self.y = x, y
        self.current_point = self
        self.next_point = None

    def _bounds(self):
        """ 
        Coordinatess of rectangle surrounding circlular object.
        """
        return (self.x + self._POINT_RADIUS*self._COS_0,   self.y + self._POINT_RADIUS*self._SIN_270,
                self.x + self._POINT_RADIUS*self._COS_180, self.y + self._POINT_RADIUS*self._SIN_90)

    def _draw(self):
        self.oval = self._graph.canvas.create_oval(self._bounds(), fill=self._color, width=0)

    def distance_to_point(self, point):
        '''
        The real distance between self and another point.
        '''
        return math.sqrt((self.x-point.x)**2 + (self.y-point.y)**2)

    def angle_to_point(self, point):
        '''
        The angle (in degrees) from self to another point, adjusted so that is a multiple of 10 (so that we don't end up with
        values that are slightly off when checking angle equality).
        '''
        if self == point: 
            return 0.0
        angle = math.ceil(180/math.pi * (math.atan2(point.y-self.y, self.x-point.x)-math.pi/2))
        if angle < 0:
            angle += math.fabs(angle) % 10
        elif angle > 0:
            angle -= math.fabs(angle) % 10
        return angle

    def __str__(self):
        return "(" + str(self.x) + ", " + str(self.y) + ")"

    def __eq__(self, other):
        return self.current_point.x == other.current_point.x and self.current_point.y == other.current_point.y

    def __ne__(self, other):
        return not self.__eq__(other)

class Vertex(Point):
    '''
    A Vertex in a graph.
    '''
    def __init__(self, graph, x, y):
        super(Vertex, self).__init__(graph, x, y)
        self._draw()

class Robot(Point):
    '''
    A robot that can move along a graph.
    '''
    _ROBOT_COLORS = ["red", "blue", "green", "cyan", "yellow", "magenta"]
    _SPEED = 5

    def __init__(self, graph, x, y, i):
        super(Robot, self).__init__(graph, x, y)
        self._color = self._ROBOT_COLORS[i%len(self._ROBOT_COLORS)]

        # A boolean value that is True if the robot's INTENDED position is the center of the graph, and False otherwise.
        # Note that this will be False if the robot is merely passing through the center to get somewhere else.
        self.centered = True

        # We make robots slightly smaller than vertices.
        self._POINT_RADIUS -= 1
        self._draw()

    def move_to(self, point):
        '''
        Move this robot to point.
        '''
        if self.next_point is None:
            self.next_point = point
        distance = self.distance_to_point(self.next_point)
        angle = self.angle_to_point(self.next_point)
        if self.current_point != point:
            self._graph.canvas.move(self.oval, self._SPEED*sin(angle), self._SPEED*cos(angle))
        x0, y0, x1, y1 = self._graph.canvas.coords(self.oval)
        self.x = (x0+x1)/2
        self.y = (y0+y1)/2
        if distance <= self._SPEED:
            self.next_point = None
            self.current_point = point
            self.x, self.y = self.current_point.x, self.current_point.y

class Edge:
    '''
    An edge in a graph, with a start and end-point.
    '''
    def __init__(self, graph, startVertex, endVertex):
        self._graph, self.start, self.end = graph, startVertex, endVertex
        self._graph.canvas.create_line(startVertex.x, startVertex.y, endVertex.x, endVertex.y)

    def midpoint(self):
        '''
        A Point that lies in the middle of an edge, so that we can move a Robot on an edge.
        '''
        return Point(self._graph, (self.start.x+self.end.x)/2, (self.start.y+self.end.y)/2)

    def __str__(self):
        return str(self.start) + " to " + str(self.end)
