from maass_sqfreelevel_computations import find_single_ev_linearized, find_evs, group_data
import sys
import time

def main():
    if len(sys.argv) < 5:
        print("  Usage:   sage progname N R   radius symtype")
        print("  example: sage progname 1 9.5 0.5 -1")
        sys.exit()
    level = int(sys.argv[1])
    R = float(sys.argv[2])
    radius = float(sys.argv[3])
    symtype = int(sys.argv[4])
    gd = group_data(level)
    print("Searching for an eigenvalue in B({}, {}) with symtype {}".format(R, 2*radius, symtype))
    val = find_single_ev_linearized(R, radius, None, gd, symtype, verbosity=2, allsigns=True)
    if val != None:
        print("{} is an eigenvalue.".format(val[0]))
    else:
        print("None found.")


def other_main():
    if len(sys.argv) < 5:
        print("  Usage:   sage progname N R1 R2 diff symtype")
        print("  example: sage progname 1 9  11 0.1  -1")
        sys.exit()
    level = int(sys.argv[1])
    R1 = float(sys.argv[2])
    R2 = float(sys.argv[3])
    diff = float(sys.argv[4])
    symtype = int(sys.argv[5])

    start = time.time()
    gd = group_data(level)
    evs = find_evs(R1, R2, None, gd, diff=diff, symmetry=symtype, verbosity=2, allsigns=True)

    print("## Finished in {} seconds".format(time.time() - start))
    if evs:
        prev = evs[0]
        print('#----- ev -----')
        print(prev[0])
        print(prev[2], "symtype {}".format(symtype))
        for val in prev[1]:
            print("  {}".format(val))
        print()

        for ev in evs[1:]:
            #if abs(ev[0] - prev[0]) < 1e-4:
            #    continue
            print('#----- ev -----')
            print(ev[0])
            print(ev[2], "symtype {}".format(symtype))
            for val in ev[1]:
                print("  {}".format(val))
            print()
            prev = ev

if __name__ == "__main__":
    other_main()
    # main()
