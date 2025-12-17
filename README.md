Starter code and data for traveling salesman problem

## Files in this directory

* datareader.cpp : example code to read in the data files (use Makefile)
* datareader.py  : example code to read in the data files
* cities23.dat   : list of coordinates for 23 cities in North America
* cities150.dat  : 150 cities in North America
* cities1k.dat   : 1207 cities in North America
* cities2k.dat   : 2063 cities around the world (PHYS 5630 only; not used for PHYS 3630)
* routeplot.py   : original code to plot the globe and salesman's path
* routeplot2.py  : updated plotter compatible with `sales2.py` outputs

Route plotting usage (recommended):
python routeplot2.py cities150.dat cities150_opt.dat --region NA

---

## Method summary

This project uses simulated annealing to minimize the total closed-loop tour length (km) for a set of cities on Earth using great-circle (haversine) distance.

Trial move method:
At each Metropolis step, two cities are selected uniformly at random and their positions in the tour are swapped to generate a trial configuration.

Acceptance rule:
Metropolis criterion: accept if ΔL ≤ 0; otherwise accept with probability exp(-ΔL / T).

Initial temperature:
For each dataset, the initial temperature is estimated automatically (`--autoT0`) from sampled uphill moves to target a high initial acceptance probability.

Performance optimization (sales2.py):
To make the 1k-city run fast, `sales2.py` uses an O(1) local ΔL update for swap moves instead of recomputing the full O(N) tour length each step.

---

## Results (PHYS 3630)

### Dataset: cities23.dat

Initial distance: 48148.58 km  
Optimized distance: 13404.55 km  
Runtime: 36.496 s  
Estimated T0 (autoT0): 9053.723  

Raw output (for reference):
- INITIAL_LENGTH 48148.583739
- OPTIMIZED_LENGTH 13404.554537
- EXECUTION_TIME 36.496442

Outputs:
- Optimized route file: cities23_opt.dat
- Annealing schedule plot: ancities23.png

Command line:
python sales2.py cities23.dat --seed 1 --autoT0 --alpha 0.995 --steps_per_T 2000 --Tmin 1e-3


### Dataset: cities150.dat

Initial distance: 333808.23 km  
Optimized distance: 55063.25 km  
Runtime: 36.745 s  
Estimated T0 (autoT0): 8577.579  

Raw output (for reference):
- INITIAL_LENGTH 333808.230020
- OPTIMIZED_LENGTH 55063.252934
- EXECUTION_TIME 36.745272

Outputs:
- Optimized route file: cities150_opt.dat
- Annealing schedule plot: an150.png

Command line:
python sales2.py cities150.dat --seed 1 --autoT0 --alpha 0.995 --steps_per_T 2000 --Tmin 1e-3


### Dataset: cities1k.dat

Initial distance: 2730084.98 km  
Optimized distance: 273222.22 km  
Runtime: 37.999 s  
Estimated T0 (autoT0): 8151.241  

Raw output (for reference):
- INITIAL_LENGTH 2730084.978602
- OPTIMIZED_LENGTH 273222.221704
- EXECUTION_TIME 37.999447

Outputs:
- Optimized route file: cities1k_opt.dat
- Annealing schedule plot: an1k.png

Command line:
python sales2.py cities1k.dat --seed 1 --autoT0 --alpha 0.995 --steps_per_T 2000 --Tmin 1e-3 --schedule_stride 10

---

## How to run

General command format:
python sales2.py <cities.dat> --seed 1 --autoT0 --alpha 0.995 --steps_per_T 2000 --Tmin 1e-3

Examples:
python sales2.py cities23.dat  --seed 1 --autoT0 --alpha 0.995 --steps_per_T 2000 --Tmin 1e-3
python sales2.py cities150.dat --seed 1 --autoT0 --alpha 0.995 --steps_per_T 2000 --Tmin 1e-3
python sales2.py cities1k.dat  --seed 1 --autoT0 --alpha 0.995 --steps_per_T 2000 --Tmin 1e-3 --schedule_stride 10

Create logs:
mkdir -p logs
for f in cities23.dat cities150.dat cities1k.dat; do
  python -u sales2.py "$f" --seed 1 --autoT0 --alpha 0.995 --steps_per_T 2000 --Tmin 1e-3 --schedule_stride 10 \
    > "logs/${f%.dat}.log" 2>&1
done




