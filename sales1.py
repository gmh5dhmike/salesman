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


def haversine(lon1_deg, lat1_deg, lon2_deg, lat2_deg):
    """Great-circle distance between two points on Earth, in km."""
    lon1 = math.radians(lon1_deg)
    lat1 = math.radians(lat1_deg)
    lon2 = math.radians(lon2_deg)
    lat2 = math.radians(lat2_deg)

    dlon = lon2 - lon1
    dlat = lat2 - lat1

    a = math.sin(dlat / 2.0) ** 2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2.0) ** 2
    c = 2.0 * math.atan2(math.sqrt(a), math.sqrt(1.0 - a))
    return R_EARTH * c


def tour_length(cities, tour):
    """
    Total loop length for a tour (including return to start).
    tour: list of indices into cities.
    """
    total = 0.0
    n = len(tour)
    for i in range(n):
        i1 = tour[i]
        i2 = tour[(i + 1) % n]  # wrap
        lon1, lat1, _ = cities[i1]
        lon2, lat2, _ = cities[i2]
        total += haversine(lon1, lat1, lon2, lat2)
    return total


def propose_swap(tour):
    """Return a new tour obtained by swapping two random positions."""
    n = len(tour)
    i, j = random.sample(range(n), 2)
    new_tour = tour[:]
    new_tour[i], new_tour[j] = new_tour[j], new_tour[i]
    return new_tour


def estimate_T0(cities, tour, samples=300, p0=0.8):
    """
    Estimate an initial temperature T0 so that a typical uphill move
    is accepted with probability ~ p0 at the start.

    We sample random swap moves from the current tour, collect positive
    cost increases dE > 0, and set:

        T0 = - mean(dE_pos) / ln(p0)

    If no uphill moves are observed (rare), fall back to a small T0.
    """
    current_L = tour_length(cities, tour)
    dE_pos = []

    for _ in range(samples):
        cand = propose_swap(tour)
        new_L = tour_length(cities, cand)
        dE = new_L - current_L
        if dE > 0:
            dE_pos.append(dE)

    if not dE_pos:
        return 1.0

    mean_dE = sum(dE_pos) / len(dE_pos)
    return -mean_dE / math.log(p0)


def metropolis_step(cities, tour, current_length, T):
    """
    One Metropolis step at temperature T.
    Returns (new_tour, new_length).
    """
    candidate = propose_swap(tour)
    new_length = tour_length(cities, candidate)
    dE = new_length - current_length

    if dE <= 0:
        return candidate, new_length
    else:
        if random.random() < math.exp(-dE / T):
            return candidate, new_length
        else:
            return tour, current_length


def simulated_annealing(cities, initial_tour,
                        T0=2000.0, alpha=0.995,
                        steps_per_T=2000, Tmin=1e-3):
    """
    Basic simulated annealing loop.
    Objective function = total path length (km).
    Returns (best_tour, best_length, schedule),
    where schedule is a list of (T, best_length) pairs.
    """
    tour = initial_tour[:]
    current_length = tour_length(cities, tour)
    best_tour = tour[:]
    best_length = current_length

    schedule = []  # (temperature, best_length)

    T = T0
    while T > Tmin:
        for _ in range(steps_per_T):
            tour, current_length = metropolis_step(cities, tour, current_length, T)
            if current_length < best_length:
                best_length = current_length
                best_tour = tour[:]
        schedule.append((T, best_length))
        T *= alpha

    return best_tour, best_length, schedule


def write_tour_lonlat(filename, cities, tour):
    """
    Write an ordered list of city coordinates for routeplot.py.
    Only longitude and latitude columns, no names.
    """
    with open(filename, "w") as f:
        f.write("# optimized tour: lon lat\n")
        for idx in tour:
            lon, lat, _ = cities[idx]
            f.write(f"{lon: .6f} {lat: .6f}\n")


def plot_schedule(input_file, schedule):
    """
    Plot annealing schedule: total distance vs temperature.
    Output filename follows an[150,1k,2k].png if we can detect the size,
    otherwise an_<basename>.png.
    """
    if not schedule:
        return

    temps = [T for (T, L) in schedule]
    dists = [L for (T, L) in schedule]

    base, _ = os.path.splitext(os.path.basename(input_file))

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
    plt.ylabel("Total distance (km)")
    plt.title(f"Annealing schedule for {input_file}")
    plt.grid(True)
    plt.savefig(outname, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Wrote annealing schedule plot to {outname}")


def main():
    print(">> sales1.py starting up")
    parser = argparse.ArgumentParser(description="Simulated annealing TSP solver")
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
    parser.add_argument("--constraints", action="store_true",
                        help="Apply additional path constraints (not used in base run)")
    parser.add_argument("--autoT0", action="store_true",
                        help="Estimate initial temperature from trial swap moves")

    args = parser.parse_args()
    random.seed(args.seed)

    cities = read_cities(args.input_file)
    N = len(cities)
    print(f"Read {N} cities from {args.input_file}")

    # Initial random tour
    initial_tour = list(range(N))
    random.shuffle(initial_tour)

    initial_length = tour_length(cities, initial_tour)
    print(f"Initial tour length: {initial_length:.2f} km")

    # Choose T0
    T0 = args.T0
    if args.autoT0:
        T0 = estimate_T0(cities, initial_tour)
        print(f"Estimated initial temperature T0 = {T0:.3f}")

    # Time the annealing
    t0 = time.time()
    best_tour, best_length, schedule = simulated_annealing(
        cities,
        initial_tour,
        T0=T0,
        alpha=args.alpha,
        steps_per_T=args.steps_per_T,
        Tmin=args.Tmin,
    )
    t1 = time.time()
    exec_time = t1 - t0

    print(f"Optimized tour length: {best_length:.2f} km")
    print(f"Execution time: {exec_time:.3f} s")

    # Write optimized path file for routeplot.py
    base, _ = os.path.splitext(args.input_file)
    optfile = base + "_opt.dat"
    write_tour_lonlat(optfile, cities, best_tour)
    print(f"Wrote optimized route to {optfile}")

    print(f"INITIAL_LENGTH {initial_length:.6f}")
    print(f"OPTIMIZED_LENGTH {best_length:.6f}")
    print(f"EXECUTION_TIME {exec_time:.6f}")

    # Annealing schedule plot
    plot_schedule(args.input_file, schedule)


if __name__ == "__main__":
    main()
