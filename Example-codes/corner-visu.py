import numpy as np
import math

X = 6
Y = 6
HALO = 2


def print_domain(array):
    np.set_printoptions(precision=3)
    print(array)


def clear_halos(array, x=X, y=Y, halo_size=HALO):
    for i in range(0, halo_size):
        for j in range(0, y + 2 * halo_size):
            array[i, j] = 0
    for i in range(x + halo_size, x + 2 * halo_size):
        for j in range(0, y + 2 * halo_size):
            array[i, j] = 0

    for j in range(0, halo_size):
        for i in range(0, x + 2 * halo_size):
            array[i, j] = 0
    for j in range(y + halo_size, y + 2 * halo_size):
        for i in range(0, x + 2 * halo_size):
            array[i, j] = 0


def fill4corneres(array, direction=1):
    xstart = HALO
    ystart = HALO
    xend = X + HALO - 1
    yend = Y + HALO - 1
    if direction == 1:
        # sw corner
        array[xstart, ystart] = array[xstart, ystart + 1]
        array[xstart - 1, ystart] = array[xstart, ystart + 2]
        # se corner
        array[xend + 1, ystart] = array[xend, ystart + 2]
        array[xend, ystart] = array[xend, ystart + 1]
        # ne corner
        array[xend, yend] = array[xend, yend - 1]
        array[xend + 1, yend] = array[xend, yend - 2]
        # nw corner
        array[xstart, yend] = array[xstart, yend - 1]
        array[xstart - 1, yend] = array[xstart, yend - 2]
    else:
        # sw corner
        array[xstart, ystart] = array[xstart + 1, ystart]
        array[xstart, ystart - 1] = array[xstart + 2, ystart]
        # se corner
        array[xend, ystart] = array[xend - 1, ystart]
        array[xend, ystart - 1] = array[xend - 2, ystart]
        # ne corner
        array[xend, yend] = array[xend - 1, yend]
        array[xend, yend + 1] = array[xend - 2, yend]
        # nw corner
        array[xstart, yend] = array[xstart + 1, yend]
        array[xstart, yend + 1] = array[xstart + 2, yend]


def great_circle_dist(q1, q2, radius=1):
    beta = (
        math.asin(
            math.sqrt(
                math.sin((q1[1] - q2[1]) / 2.0) ** 2
                + math.cos(q1[1]) * math.cos(q2[1]) * math.sin((q1[0] - q2[0]) / 2.0) ** 2
            )
        )
        * 2.0
    )
    return beta * radius


def extrap_corner(p0, p1, p2, q1, q2):
    x1 = great_circle_dist(p1, p0)
    x2 = great_circle_dist(p2, p0)

    return q1 + x1 / (x2 - x1) * (q1 - q2)


def a2b_ord4(in_field, out_field):
    xstart = HALO
    ystart = HALO
    xend = X + HALO - 1
    yend = Y + HALO - 1
    # sw corner
    tmp = in_field
    out_field[xstart + 1, ystart + 1] = tmp


if __name__ == "__main__":
    print("Visualization of corner handeling\n")

    y = np.arange((X + 2 * HALO) * (Y + 2 * HALO)).reshape(X + 2 * HALO, Y + 2 * HALO)
    print_domain(y)

    clear_halos(y)
    print_domain(y)

    fill4corneres(y)
    fill4corneres(y, 2)
    print_domain(y)
