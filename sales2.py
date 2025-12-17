#!/usr/bin/env python3
import math
import random
import argparse
import os
import time
import matplotlib.pyplot as plt

R_EARTH = 6371.0  # km


def read_cities(filename):
    """
    Read cities from a data file.
    Expected line format (ignoring comments starting with '#'):
        longitude latitude "City Name"
    Returns list of (lon_deg, lat_deg, name).
    """
    cities = []
    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split(maxsplit=2)
            lon = float(parts[0])
            lat = float(parts[1])
            name = parts[2].strip().strip('"') if len(parts) > 2 else ""
            cities.append((lon, lat, name))
    return cities


def preprocess_cities(cities):
    """
    Precompute lon/lat in radians and cos(lat) to speed up distance calculations.
    Returns list of tuples: (lon_rad, lat_rad, cos_lat, name)
    """
    out = []
    for lon_deg, lat_deg, name in cities:
        lon = math.radians(lon_deg)
        lat = math.radians(lat_deg)
        out.append((lon, lat, math.cos(lat), name))
    return out


def haversine_rad(c1, c2):
    """
    Great-circle distance (km) between two preprocessed city records:
      c = (lon_rad, lat_rad, cos_lat, name)
    """
    lon1, lat1, cos1, _ = c1
    lon2, lat2, cos2, _ = c2

    dlon = lon2 - lon1
    dlat = lat2 - lat1

    sin_dlat = math.sin(dlat * 0.5)
    sin_dlon = math.sin(dlon * 0.5)

    a = sin_dlat * sin_dlat + cos1 * cos2 * (sin_dlon * sin_dlon)
    c = 2.0 * math.atan2(math.sqrt(a), math.sqrt(1.0 - a))
    return R_EARTH * c


def edge_len(pcities, a_idx, b_idx):
    """Distance (km) between two cities indexed by a_idx, b_idx using preprocessed cities."""
    return haversine_rad(pcities[a_idx], pcities[b_idx])


def tour_length(pcities, tour):
    """
    Total loop length for a tour (including return to start).
    tour: list of indices into pcities.
    """
    total = 0.0
    n = len(tour)
    for k in range(n):
        a = tour[k]
        b = tour[(k + 1) % n]
        total += edge_len(pcities, a, b)
    return total


def estimate_T0(pcities, tour, samples=300, p0=0.8):
    """
    Estimate an initial temperature T0 so that a typical uphill move
    is accepted with probability ~ p0 at the start.

    T0 = - mean(dE_pos) / ln(p0)
    where dE_pos are sampled positive cost increases.
    """
    current_L = tour_length(pcities, tour)
    dE_pos = []
    n = len(tour)

    for _ in range(samples):
        i, j = random.sample(range(n), 2)
        dE = delta_len_swap(pcities, tour, i, j)
        if dE > 0:
            dE_pos.append(dE)

    if not dE_pos:
        return 1.0

    mean_dE = sum(dE_pos) / len(dE_pos)
    return -mean_dE / math.log(p0)


def delta_len_swap(pcities, tour, i, j):
    """
    O(1) change in tour length if we swap the *positions* i and j in the tour.
    Handles adjacency and wrap-around correctly.
    """
    n = len(tour)
    if i == j:
        return 0.0
    if i > j:
        i, j = j, i

    a = tour[i]
    b = tour[j]

    im1 = tour[(i - 1) % n]
    ip1 = tour[(i + 1) % n]
    jm1 = tour[(j - 1) % n]
    jp1 = tour[(j + 1) % n]

    # Adjacent in tour order (i directly before j)
    if (i + 1) % n == j:
        # ... im1 - a - b - jp1 ...
        old = edge_len(pcities, im1, a) + edge_len(pcities, a, b) + edge_len(pcities, b, jp1)
        new = edge_len(pcities, im1, b) + edge_len(pcities, b, a) + edge_len(pcities, a, jp1)
        return new - old

    # Also handle the wrap adjacency case where i=0, j=n-1 (they are adjacent via wrap)
    if i == 0 and j == n - 1:
        # ... jm1 - b - a - ip1 ... (since b is last, a is first)
        old = edge_len(pcities, jm1, b) + edge_len(pcities, b, a) + edge_len(pcities, a, ip1)
        new = edge_len(pcities, jm1, a) + edge_len(pcities, a, b) + edge_len(pcities, b, ip1)
        return new - old

    # Non-adjacent case: remove 4 edges and add 4 edges
    old = (
        edge_len(pcities, im1, a) + edge_len(pcities, a, ip1) +
        edge_len(pcities, jm1, b) + edge_len(pcities, b, jp1)
    )
    new = (
        edge_len(pcities, im1, b) + edge_len(pcities, b, ip1) +
        edge_len(pcities, jm1, a) + edge_len(pcities, a, jp1)
    )
    return new - old


def metropolis_step_inplace(pcities, tour, current_length, T):
    """
    One Metropolis step at temperature T using:
      - in-place swap (no tour copy)
      - O(1) delta length update (no full tour recompute)
    Returns (tour, new_length).
    """
    n = len(tour)
    i, j = random.sample(range(n), 2)

    dE = delta_len_swap(pcities, tour, i, j)

    # Metropolis acceptance
    if dE <= 0 or random.random() < math.exp(-dE / T):
        tour[i], tour[j] = tour[j], tour[i]
        return tour, current_length + dE
    else:
        return tour, current_length


def simulated_annealing(pcities, initial_tour,
                        T0=2000.0, alpha=0.995,
                        steps_per_T=2000, Tmin=1e-3,
                        schedule_stride=1):
    """
    Simulated annealing loop.

    schedule_stride controls how often we record the schedule point:
      - 1 records every temperature
      - 10 records every 10 temperatures (faster plotting / smaller list)
    """
    tour = initial_tour[:]  # keep a working copy
    current_length = tour_length(pcities, tour)
    best_tour = tour[:]
    best_length = current_length

    schedule = []
    T = T0
    tcount = 0

    while T > Tmin:
        for _ in range(steps_per_T):
            tour, current_length = metropolis_step_inplace(pcities, tour, current_length, T)
            if current_length < best_length:
                best_length = current_length
                best_tour = tour[:]

        if schedule_stride <= 1 or (tcount % schedule_stride == 0):
            schedule.append((T, best_length))

        T *= alpha
        tcount += 1

    return best_tour, best_length, schedule


def write_tour_lonlat(filename, cities_deg, tour):
    """
    Write an ordered list of city coordinates for routeplot.py.
    Only longitude and latitude columns, no names.
    """
    with open(filename, "w") as f:
        f.write("# optimized tour: lon lat\n")
        for idx in tour:
            lon, lat, _ = cities_deg[idx]
            f.write(f"{lon:.6f} {lat:.6f}\n")


def plot_schedule(input_file, schedule):
    """
    Plot annealing schedule: best distance vs temperature.
    """
    if not schedule:
        return

    temps = [T for (T, L) in schedule]
    dists = [L for (T, L) in schedule]

    base = os.path.splitext(os.path.basename(input_file))[0]

    if "150" in base:
        tag = "150"
    elif "1k" in base:
        tag = "1k"
    elif "2k" in base:
        tag = "2k"
    else:
        tag = base

    outname = f"an{tag}.png"

    plt.figure()
    plt.plot(temps, dists, marker="o", markersize=3)
    plt.xlabel("Temperature")
    plt.ylabel("Best distance (km)")
    plt.title(f"Annealing schedule for {input_file}")
    plt.grid(True)
    plt.savefig(outname, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Wrote annealing schedule plot to {outname}")


def main():
    print(">> sales2.py starting up")
    parser = argparse.ArgumentParser(description="Simulated annealing TSP solver (optimized swap delta)")
    parser.add_argument("input_file", help="City data file, e.g. cities23.dat")
    parser.add_argument("--T0", type=float, default=2000.0,
                        help="Initial temperature (ignored if --autoT0 is set)")
    parser.add_argument("--alpha", type=float, default=0.995,
                        help="Cooling factor (0 < alpha < 1)")
    parser.add_argument("--steps_per_T", type=int, default=2000,
                        help="Metropolis steps per temperature")
    parser.add_argument("--Tmin", type=float, default=1e-3,
                        help="Final temperature cutoff")
    parser.add_argument("--seed", type=int, default=1,
                        help="Random seed")
    parser.add_argument("--autoT0", action="store_true",
                        help="Estimate initial temperature from trial swap moves")
    parser.add_argument("--schedule_stride", type=int, default=1,
                        help="Record schedule every k temperatures (default 1)")
    args = parser.parse_args()

    random.seed(args.seed)

    cities_deg = read_cities(args.input_file)
    pcities = preprocess_cities(cities_deg)
    N = len(pcities)
    print(f"Read {N} cities from {args.input_file}")

    # Initial random tour
    tour = list(range(N))
    random.shuffle(tour)

    initial_length = tour_length(pcities, tour)
    print(f"Initial tour length: {initial_length:.2f} km")

    # Choose T0
    T0 = args.T0
    if args.autoT0:
        T0 = estimate_T0(pcities, tour)
        print(f"Estimated initial temperature T0 = {T0:.3f}")

    # Time the annealing
    t_start = time.time()
    best_tour, best_length, schedule = simulated_annealing(
        pcities,
        tour,
        T0=T0,
        alpha=args.alpha,
        steps_per_T=args.steps_per_T,
        Tmin=args.Tmin,
        schedule_stride=args.schedule_stride,
    )
    exec_time = time.time() - t_start

    print(f"Optimized tour length: {best_length:.2f} km")
    print(f"Execution time: {exec_time:.3f} s")

    # Write optimized path file for routeplot.py (uses original deg coords)
    base, _ = os.path.splitext(args.input_file)
    optfile = base + "_opt.dat"
    write_tour_lonlat(optfile, cities_deg, best_tour)
    print(f"Wrote optimized route to {optfile}")

    print(f"INITIAL_LENGTH {initial_length:.6f}")
    print(f"OPTIMIZED_LENGTH {best_length:.6f}")
    print(f"EXECUTION_TIME {exec_time:.6f}")

    plot_schedule(args.input_file, schedule)


if __name__ == "__main__":
    main()
