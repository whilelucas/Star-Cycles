#!/usr/bin/env python2

from sage.all import *
import cubical_complex
import itertools
import logging
import collections

# Ensure we've got chomp
from sage.interfaces.chomp import have_chomp
#assert have_chomp() is True

# Logging configuration: by default, produce no output
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# Return lookup_Y and generate_tree_Y in lookup and generate_tree respectively to get
# the Abrams-discretized model D_n(Y).
def lookup_Y(n):
    """\
    Lookup returns a list of points in R^3, it amounts to a labeling of the
    vertices of the Y graph. Examples, from left to right: ``lookup(1)``,
    ``lookup(2)``, and ``lookup(3)``:

     |         2              4
     0__       |              |
    /          1__0           3
              /               |
             3                2__1__0
                             /
                            5
                           /
                          6

    ``lookup(2)`` would be represented as a list like the following:

        [ (1,0,0) # index/vertex 0 represents the point (1,0,0)
        , (0,0,0) # index/vertex 1 represents the origin
        , (0,1,0) # index/vertex 2 ...
        , (0,0,1) # index/vertex 3 ...
        ]

    Examples:

        >>> lookup(1)
        [(0, 0, 0)]

        >>> lookup(2)
        [(0, 1, 0), (0, 0, 0), (1, 0, 0), (0, 0, 1)]

        >>> import random
        >>> n = random.randint(1, 100)
        >>> lookup_ = lookup(n)
        >>> lookup_[n-1]
        (0, 0, 0)

      * The label-to-coordinate mapping should be unique:

        >>> assert len(frozenset(lookup_)) == len(lookup_)
    """

    assert n > 0

    lookup = [(0, n - t - 1, 0) for t in range(n - 1)]
    lookup.append((0, 0, 0))
    lookup.extend([(t, 0, 0) for t in range(1, n)])
    lookup.extend([(0, 0, t) for t in range(1, n)])

    assert lookup != []

    return lookup

def generate_tree_Y(n):
    """\
    Returns a list of possible "downstream" moves (moves from one point to a
    point with a greater label). For instance, in the case of C^n(G)
    (e.g. C^2(G)),

     * Initial: if a point is at 0, it can move to 1. If a point is at 1, it
       can move to 2, ... , n-1 (e.g. 2).
     * Center: If a point is at n-1 (e.g. 2), it can move to either n or 2n-1
       (e.g. 3 or 5).
     * Legs: If a point is at n, it can move to n+1, etc. If a point is at
       2n-1, it can move to 2n, ... , 3n-3.

    Examples:

        >>> generate_tree(1) # TODO: this seems incorrect
        ([1, 1], [], [])

        >>> # The vertex at 0 can go to 1, the vertex at 1 can go to 2 or 3...
        >>> generate_tree(2)
        ([1], [2, 3], [], [])

        >>> import random
        >>> n = random.randint(1, 100)
        >>> assert generate_tree(n)[n-1] == [n, 2*n-1]
        >>> # TODO: this fails on input "1"
        >>> generate_tree(n)[2*n-2]  # end of leg 1
        []
        >>> generate_tree(n)[3*n-3]  # end of leg 2
        []
    """
    assert n > 0

    tree = []
    for point in range(n - 1):
        tree.append([point + 1])
    tree.append([n, 2 * n - 1])
    for point in range(n, 2 * n - 2):
        tree.append([point + 1])
    tree.append([])
    for point in range(2 * n - 1, 3 * n - 3):
        tree.append([point + 1])
    tree.append([])

    # Points can only move downstream, to points of greater labels
    for index, lst in enumerate(tree):
        assert filter(lambda x: x <= index, lst) == []

    return tuple(tree)

# Return lookup_X and generate_tree_X in lookup and generate_tree respectively to get
# the Abrams-discretized model D_n(X).
def lookup_X(n):
    assert n > 0

    lookup = [(0, n-t-1, 0, 0) for t in range(n-1)]
    lookup.append((0, 0, 0, 0))
    lookup.extend([(t, 0, 0, 0) for t in range(1 ,n)])
    lookup.extend([(0, 0, t, 0) for t in range(1 ,n)])
    lookup.extend([(0, 0, 0, t) for t in range (1, n)])

    assert lookup != []

    return lookup

def generate_tree_X(n):
    tree = []
    for point in range (n-1):
        tree.append([point+1])
    tree.append([n, 2*n-1, 3*n-2])
    for point in range (n , 2*n-2):
        tree.append([point+1])
    tree.append([])
    for point in range (2*n-1, 3*n-3):
        tree.append([point+1])
    tree.append([])
    for point in range(3*n-2, 4*n-4):
        tree.append([point+1])
    tree.append([])
    return tuple(tree)

# Choose either lookup_Y and 
def lookup(n):
    return lookup_Y(n)

def generate_tree(n):
    return generate_tree_Y(n) 

def iterate_over_conf(T, n):
    """\
    Enumerate all the possible configurations of points at vertices, which will
    be the 0-cells. Gives lists of length n which consist of the positions
    (vertex labels) of different points. There is probably a more efficient way
    to do this.

    Examples:
       # >>> gen = iterate_over_conf(generate_tree(2), 2)
       # >>> list(gen)[:3] # First three items
       # [[0, 1], [0, 2], [0, 3]]
    """
    # Labels table shouldn't be []
    assert T != []

    for point_config in itertools.product(xrange(len(T)), repeat=n):
        if len(set(point_config)) == n:  # points all distinct
            yield point_config


def downstream_moves(point_config, T):
    """\
    Test which points in the configuration can move. We only move 'downstream'
    (i.e. to a vertex with a greater label) to avoid duplicates. Outputs a list
    of the available moves.

    Examples:
     * If we have

                         2                                     2
                         |                                     |
       T := tree(2) =    1__0  and point_config := [0, 1] =    (b)__(a)
                        /                                     /
                       3                                     3

       then point (a) is stuck, but point (b) can move to 2 or 3. Thus:

        >>> downstream_moves([0, 1], generate_tree(2))
        [[None], [2, 3]]

     * If we have

                         2                                    (a)
                         |                                     |
       T := tree(2) =    1__0  and point_config := [2, 3] =    1__0
                        /                                     /
                       3                                    (b)

       then neither can move. Thus:

        >>> downstream_moves([2, 3], generate_tree(2))
        [[None], [None]]

       In general, points at the ends of legs can't move anywhere:

        >>> import random
        >>> n = random.randint(3, 100)
        >>> moves = downstream_moves(range(n-1) + [2*n-2], generate_tree(n))
        >>> moves[n-1]  # end of leg 1
        [None]
        >>> moves = downstream_moves(range(n-1) + [3*n-3], generate_tree(n))
        >>> moves[n-1]  # end of leg 2
        [None]

    """
    assert T != []
    assert point_config != []
    # no two points can be at the same location (in the non-2-equal setting)
    assert len(point_config) == len(set(point_config))

    output = []
    point_config_set = frozenset(point_config)  # sets have O(1) lookup (hash)
    for p in point_config:
        # We can move anywhere that there isn't already a point
        downstream = [pos for pos in T[p] if pos not in point_config_set]
        if downstream == []:
            output.append([None])
        else:
            output.append(downstream)

    assert output != []

    return output

# We "tag" the downstream cubes with their origin for debugging. If we know
# what point_config and move inspired the addition of a cube to the Abrams
# discretized configuration space, then we can more easily figure out where
# things went wrong.
MoveCube = collections.namedtuple("MoveCube", ["point_config", "move", "cube"])


def downstream_cubes(point_config, T):
    """\
    Builds the highest-dimensional cubes (in the Abrams-discretized
    configuration space) that result from performing moves at the same time.

    Examples:
     * If a point at 0 can move to 1 and a point at 2 can move to 3 or 5, then
       this adds the two 2-cubes that correspond to moving 0 to 1 and 2 to 3,
       and moving 0 to 1 and 2 to 5.

       TODO: why this output?
       TODO: this _doesn't_ only construct maximal cubes

        # >>> downstream_cubes([0, 2], generate_tree(2))
        # [[[0, 0], [0, 1], [0, 0], [1, 1], [0, 0], [0, 0]]]

    See also:
      * http://doc.sagemath.org/html/en/reference/homology/sage/homology/cubical_complex.html
    """
    assert point_config != []
    assert T != []

    lookup_ = lookup(len(point_config))

    cubes = []
    # This product contains all possible combinations of moves
    for move in itertools.product(*downstream_moves(point_config, T)):
        new_cube = []
        for (next_pos, current_pos) in zip(move, point_config):
            # This point's current position in R^3
            embedded_coords = lookup_[current_pos]
            # This point doesn't move, so the intervals are trivial

            # TODO: what's going on here??
            if next_pos is None:
                intervals = [[u, u] for u in embedded_coords]
            else:
                # This point's position in R^3 after the move
                next_embedded_coords = lookup_[next_pos]
                intervals = [
                    sorted((u, v))
                    for (u, v) in zip(embedded_coords, next_embedded_coords)
                ]
            new_cube.extend(intervals)  # (each cube has 3n intervals)

        # Make a new "tagged" cube
        cubes.append(
            MoveCube(point_config, move, cubical_complex.Cube(new_cube)))

    assert cubes != []

    return cubes


def the_complex(n, maximality_check=True, logger=logger):
    """ Build the cubical complex that is the Abrams-discretized configuration
    space of n vertices on the Y graph.

    The maximality_check _should be_ set to False by default because we
    intentionally only add maximal cells to the complex with downstream_moves.
    This is not currently the case.

    Examples:

        # TODO: this fails:
        # >>> the_complex(1)

        >>> the_complex(2)
        Cubical complex with 12 vertices and 24 cubes
        >>> the_complex(2).homology()
        {0: 0, 1: Z}

        >>> len(lookup(3))
        7
        >>> 7 * 6 * 5
        210
        >>> the_complex(3)
        Cubical complex with 210 vertices and 756 cubes
        >>> the_complex(3).homology()
        {0: 0, 1: Z^13, 2: 0, 3: 0}

    """
    assert n > 0

    T = generate_tree(n)
    cubes = []

    # If any of the points are at the ends of the legs, then the
    # generated cube will be a face of one already generated.
    for point_config in iterate_over_conf(T, n):
        logger.debug("Generating downstream_cubes for {}".format(point_config))
        downstream = downstream_cubes(point_config, T)

        # for down in downstream:
        #     for cube in cubes:
        #         if down.cube.is_face(cube.cube):
        #             try:
        #                 assert 2 * n - 2 in point_config or 3 * n - 3 in point_config
        #             except AssertionError:
        #                 print("ERRR: the first contains the second")
        #                 print("cube.point_config: {}".format(cube.point_config))
        #                 print("cube.move: {}".format(cube.move))
        #                 print("cube.cube: {}".format(cube.cube))
        #                 print("down.point_config: {}".format(down.point_config))
        #                 print("down.move: {}".format(down.move))
        #                 print("down.cube: {}".format(down.cube))
        #                 raise

        cubes.extend(downstream)

    cubes = map(lambda t: t.cube, cubes)
    downstream = map(lambda x: x.cube, downstream)

    return cubical_complex.CubicalComplex(
        cubes, maximality_check=maximality_check)
