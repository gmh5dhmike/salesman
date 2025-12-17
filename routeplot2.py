import matplotlib.pyplot as plt
import numpy as np
import argparse
import os


def load_xy(filename):
    """Load only first two columns (lon lat), skipping lines beginning with '#'."""
    data = []
    with open(filename) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.split()
            data.append([float(cols[0]), float(cols[1])])
    # close the path
    if data:
        data.append(data[0])
    return np.array(data)


def load_polygons(filename):
    """Load longitude/latitude polygons separated by blank lines."""
    polygons = []
    current = []

    with open(filename) as f:
        for line in f:
            line = line.strip()

            if not line:
                if current:
                    polygons.append(np.array(current))
                    current = []
                continue

            parts = line.split()
            lon = float(parts[0])
            lat = float(parts[1])
            current.append([lon, lat])

    if current:
        polygons.append(np.array(current))

    return polygons


def infer_region_from_name(infile, default="NA"):
    base = os.path.basename(infile)
    if "2k" in base or "world" in base.lower():
        return "World"
    return default


def make_plot(infile, optfile=None, region="NA", outpdf=None, mapfile="world_50m.dat"):
    """
    infile: required city list (original order)
    optfile: optional optimized route file produced by sales2.py (e.g., cities150_opt.dat)
    region: "NA" or "World"
    outpdf: output pdf filename (optional)
    mapfile: polygon coastline file (optional; used if present)
    """

    cities_orig = load_xy(infile)
    cities_opt = load_xy(optfile) if optfile else None

    fig, ax = plt.subplots(figsize=(10, 6))

    # Plot world outline if file exists
    if os.path.exists(mapfile):
        polygons = load_polygons(mapfile)
        for poly in polygons:
            ax.plot(poly[:, 0], poly[:, 1], color="black", lw=0.8)
    else:
        print(f"Note: map file '{mapfile}' not found. Plotting routes without coastline.")

    # Original order (thin line)
    if len(cities_orig) > 0:
        ax.plot(cities_orig[:, 0], cities_orig[:, 1], lw=1, color="red", alpha=0.5, label="original")

    # Optimized route (thicker line + points)
    if cities_opt is not None and len(cities_opt) > 0:
        ax.plot(cities_opt[:, 0], cities_opt[:, 1],
                lw=2, color="blue", marker='o', markersize=2, label="optimized")

    ax.set_title("Traveling Salesman Route")
    ax.set_xlabel("longitude")
    ax.set_ylabel("latitude")

    # Axis ranges
    if region == "World":
        ax.set_xlim(-180, 180)
        ax.set_ylim(-90, 90)
    else:
        ax.set_xlim(-180, -60)
        ax.set_ylim(10, 75)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.legend(loc="best")

    # Output filename
    if outpdf is None:
        base = os.path.splitext(os.path.basename(infile))[0]
        if optfile:
            outpdf = f"{base}_vs_{os.path.splitext(os.path.basename(optfile))[0]}.pdf"
        else:
            outpdf = f"{base}.pdf"

    plt.savefig(outpdf, format="pdf", facecolor="white", bbox_inches="tight")
    print(f"Wrote plot to {outpdf}")

    print('Close plot window (or Ctrl-C) to exit')
    try:
        plt.show()
    except KeyboardInterrupt:
        print("Interrupted with Ctrl-C, closing plot and exiting...")
        plt.close('all')


def main():
    parser = argparse.ArgumentParser(description="Route plotter for sales2.py outputs")
    parser.add_argument("infile", help="Original city file, e.g. cities150.dat")
    parser.add_argument("optfile", nargs="?", default=None,
                        help="Optimized route file, e.g. cities150_opt.dat (optional)")
    parser.add_argument("--region", choices=["NA", "World"], default=None,
                        help="Plot region (default: NA; auto-detect if not given)")
    parser.add_argument("-o", "--out", default=None,
                        help="Output PDF filename (optional)")
    parser.add_argument("--mapfile", default="world_50m.dat",
                        help="Polygon coastline data file (default: world_50m.dat)")

    args = parser.parse_args()

    region = args.region
    if region is None:
        region = infer_region_from_name(args.infile, default="NA")

    make_plot(args.infile, args.optfile, region=region, outpdf=args.out, mapfile=args.mapfile)


if __name__ == "__main__":
    main()
