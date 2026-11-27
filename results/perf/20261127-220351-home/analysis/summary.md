# Performance results: `20261127-213925-home`, `20261127-220351-home`

267 successful runs.

## 1. Single-VM baseline (compute rate)

Serial time per step and per node-update; median of the repetitions, spread in brackets.

| N | heat-ghadi: per step | heat-ghadi: ns/node | heat-lynn: per step | heat-lynn: ns/node | heat-rim: per step | heat-rim: ns/node |
|---|---|---|---|---|---|---|
| 512 | 203.0 µs (135.8 µs–1.32 ms) | 0.78 | 341.6 µs (279.9 µs–369.5 µs) | 1.31 | 139.1 µs (136.3 µs–143.0 µs) | 0.53 |
| 1024 | 599.5 µs (566.9 µs–4.17 ms) | 0.57 | 1.69 ms (1.57 ms–1.74 ms) | 1.61 | 950.1 µs (883.3 µs–951.9 µs) | 0.91 |
| 2048 | 3.84 ms (3.64 ms–17.08 ms) | 0.92 | 7.43 ms (7.29 ms–7.93 ms) | 1.77 | 4.39 ms (4.26 ms–4.40 ms) | 1.05 |
| 4096 | 15.84 ms (15.71 ms–73.87 ms) | 0.94 | 28.57 ms (27.89 ms–30.20 ms) | 1.70 | 19.24 ms (19.15 ms–19.25 ms) | 1.15 |
| 8192 | 63.88 ms (63.47 ms–181.13 ms) | 0.95 | 116.90 ms (115.53 ms–118.52 ms) | 1.74 | 79.75 ms (79.11 ms–81.77 ms) | 1.19 |

## Model parameters

Network from ping-pong folder `results/pingpong/20261127-205307`; model mode `pair`.

| VM pair | t_s | t_w | bandwidth |
|---|---|---|---|
| heat-ghadi <-> heat-rim | 127.86 ms | 3.030e-06 s/B | 0.330 MB/s |
| heat-lynn <-> heat-ghadi | 22.28 ms | 9.981e-07 s/B | 1.002 MB/s |
| heat-lynn <-> heat-rim | 59.64 ms | 1.404e-05 s/B | 0.071 MB/s |

* t_s, t_w within a VM: 0.3 µs, 1.117e-10 s/B
* t_node of each VM: median serial time per node per step from the baseline (table 1)

## 2. Strong scaling

Time per step (median, spread min–max), speed-up S = T1/TP, efficiency E = S/P, communication share t_comm/t_total. T1 = serial run on the first host.

### N = 512  (T1 = 417.8 µs per step on heat-lynn; fastest VM heat-rim: 139.1 µs per step)

| decomp | map | P | time/step (spread) | S | E | S_fast | comm share | model time/step | model S |
|---|---|---|---|---|---|---|---|---|---|
| serial | - | 1 | 417.8 µs (145.1 µs–630.8 µs) | 1.00 | 1.00 | 0.33 | 0% | 341.6 µs | 1.00 |
| blocks | node | 2 | 64.10 ms (63.47 ms–70.17 ms) | 0.01 | 0.00 | 0.00 | 100% | 52.89 ms | 0.01 |
| blocks | node | 3 | 98.39 ms (97.03 ms–124.56 ms) | 0.00 | 0.00 | 0.00 | 100% | 280.60 ms | 0.00 |
| blocks | node | 4 | 153.12 ms (152.41 ms–235.11 ms) | 0.00 | 0.00 | 0.00 | 100% | 353.46 ms | 0.00 |
| blocks | node | 5 | 125.98 ms (119.94 ms–187.87 ms) | 0.00 | 0.00 | 0.00 | 100% | 280.56 ms | 0.00 |
| blocks | node | 6 | 186.58 ms (166.54 ms–233.98 ms) | 0.00 | 0.00 | 0.00 | 100% | 532.14 ms | 0.00 |
| blocks | slot | 2 | 240.6 µs (232.9 µs–379.4 µs) | 1.74 | 0.87 | 0.58 | 28% | 172.6 µs | 1.98 |
| blocks | slot | 3 | 61.65 ms (55.36 ms–63.07 ms) | 0.01 | 0.00 | 0.00 | 100% | 52.84 ms | 0.01 |
| blocks | slot | 4 | 61.46 ms (55.72 ms–67.21 ms) | 0.01 | 0.00 | 0.00 | 100% | 48.74 ms | 0.01 |
| blocks | slot | 5 | 89.83 ms (83.23 ms–155.47 ms) | 0.00 | 0.00 | 0.00 | 100% | 280.56 ms | 0.00 |
| blocks | slot | 6 | 88.83 ms (64.63 ms–149.43 ms) | 0.00 | 0.00 | 0.00 | 100% | 268.18 ms | 0.00 |
| strips | node | 2 | 70.83 ms (61.10 ms–214.31 ms) | 0.01 | 0.00 | 0.00 | 100% | 52.93 ms | 0.01 |
| strips | node | 3 | 103.65 ms (60.62 ms–154.38 ms) | 0.00 | 0.00 | 0.00 | 100% | 280.70 ms | 0.00 |
| strips | node | 4 | 95.49 ms (95.40 ms–125.14 ms) | 0.00 | 0.00 | 0.00 | 100% | 280.67 ms | 0.00 |
| strips | node | 5 | 126.32 ms (104.43 ms–208.97 ms) | 0.00 | 0.00 | 0.00 | 100% | 280.65 ms | 0.00 |
| strips | node | 6 | 125.63 ms (124.02 ms–204.36 ms) | 0.00 | 0.00 | 0.00 | 100% | 280.64 ms | 0.00 |
| strips | slot | 2 | 373.0 µs (311.8 µs–860.5 µs) | 1.12 | 0.56 | 0.37 | 46% | 172.6 µs | 1.98 |
| strips | slot | 3 | 55.49 ms (54.81 ms–70.06 ms) | 0.01 | 0.00 | 0.00 | 100% | 52.87 ms | 0.01 |
| strips | slot | 4 | 58.97 ms (54.65 ms–84.28 ms) | 0.01 | 0.00 | 0.00 | 100% | 52.84 ms | 0.01 |
| strips | slot | 5 | 77.47 ms (63.99 ms–115.31 ms) | 0.01 | 0.00 | 0.00 | 100% | 280.65 ms | 0.00 |
| strips | slot | 6 | 69.41 ms (60.71 ms–86.24 ms) | 0.01 | 0.00 | 0.00 | 100% | 280.64 ms | 0.00 |

### N = 2048  (T1 = 6.65 ms per step on heat-lynn; fastest VM heat-ghadi: 3.84 ms per step)

| decomp | map | P | time/step (spread) | S | E | S_fast | comm share | model time/step | model S |
|---|---|---|---|---|---|---|---|---|---|
| serial | - | 1 | 6.65 ms (5.64 ms–7.59 ms) | 1.00 | 1.00 | 0.58 | 0% | 7.43 ms | 1.00 |
| blocks | node | 2 | 103.88 ms (86.62 ms–113.94 ms) | 0.06 | 0.03 | 0.04 | 96% | 80.97 ms | 0.09 |
| blocks | node | 3 | 113.09 ms (77.07 ms–119.69 ms) | 0.06 | 0.02 | 0.03 | 99% | 357.42 ms | 0.02 |
| blocks | node | 4 | 212.46 ms (191.98 ms–214.13 ms) | 0.03 | 0.01 | 0.02 | 99% | 700.37 ms | 0.01 |
| blocks | node | 5 | 158.52 ms (148.30 ms–161.61 ms) | 0.04 | 0.01 | 0.02 | 99% | 580.72 ms | 0.01 |
| blocks | node | 6 | 207.58 ms (178.84 ms–238.44 ms) | 0.03 | 0.01 | 0.02 | 99% | 639.38 ms | 0.01 |
| blocks | slot | 2 | 4.83 ms (4.78 ms–4.90 ms) | 1.38 | 0.69 | 0.79 | 17% | 3.72 ms | 2.00 |
| blocks | slot | 3 | 95.67 ms (85.77 ms–115.20 ms) | 0.07 | 0.02 | 0.04 | 97% | 79.73 ms | 0.09 |
| blocks | slot | 4 | 75.31 ms (73.73 ms–96.77 ms) | 0.09 | 0.02 | 0.05 | 97% | 62.78 ms | 0.12 |
| blocks | slot | 5 | 95.97 ms (78.09 ms–105.89 ms) | 0.07 | 0.01 | 0.04 | 98% | 356.43 ms | 0.02 |
| blocks | slot | 6 | 83.69 ms (72.49 ms–105.81 ms) | 0.08 | 0.01 | 0.05 | 98% | 306.59 ms | 0.02 |
| strips | node | 2 | 94.47 ms (79.75 ms–94.94 ms) | 0.07 | 0.04 | 0.04 | 96% | 81.00 ms | 0.09 |
| strips | node | 3 | 111.28 ms (96.94 ms–163.94 ms) | 0.06 | 0.02 | 0.03 | 99% | 357.51 ms | 0.02 |
| strips | node | 4 | 154.50 ms (143.96 ms–297.46 ms) | 0.04 | 0.01 | 0.02 | 99% | 581.54 ms | 0.01 |
| strips | node | 5 | 168.45 ms (163.31 ms–411.25 ms) | 0.04 | 0.01 | 0.02 | 99% | 581.17 ms | 0.01 |
| strips | node | 6 | 161.64 ms (150.99 ms–204.08 ms) | 0.04 | 0.01 | 0.02 | 99% | 580.92 ms | 0.01 |
| strips | slot | 2 | 4.49 ms (4.42 ms–5.20 ms) | 1.48 | 0.74 | 0.86 | 10% | 3.72 ms | 2.00 |
| strips | slot | 3 | 86.19 ms (76.85 ms–87.30 ms) | 0.08 | 0.03 | 0.04 | 97% | 79.76 ms | 0.09 |
| strips | slot | 4 | 76.51 ms (71.85 ms–86.29 ms) | 0.09 | 0.02 | 0.05 | 97% | 79.14 ms | 0.09 |
| strips | slot | 5 | 95.18 ms (76.22 ms–124.69 ms) | 0.07 | 0.01 | 0.04 | 99% | 356.52 ms | 0.02 |
| strips | slot | 6 | 97.00 ms (82.92 ms–223.28 ms) | 0.07 | 0.01 | 0.04 | 99% | 356.28 ms | 0.02 |

### N = 8192  (T1 = 125.29 ms per step on heat-lynn; fastest VM heat-ghadi: 63.88 ms per step)

| decomp | map | P | time/step (spread) | S | E | S_fast | comm share | model time/step | model S |
|---|---|---|---|---|---|---|---|---|---|
| serial | - | 1 | 125.29 ms (114.57 ms–126.18 ms) | 1.00 | 1.00 | 0.51 | 0% | 116.90 ms | 1.00 |
| blocks | node | 2 | 359.47 ms (323.82 ms–369.02 ms) | 0.35 | 0.17 | 0.18 | 91% | 233.83 ms | 0.50 |
| blocks | node | 3 | 341.52 ms (304.95 ms–405.54 ms) | 0.37 | 0.12 | 0.19 | 94% | 691.73 ms | 0.17 |
| blocks | node | 4 | 490.95 ms (488.41 ms–681.79 ms) | 0.26 | 0.06 | 0.13 | 96% | 2.108 s | 0.06 |
| blocks | node | 5 | 715.34 ms (640.12 ms–742.29 ms) | 0.18 | 0.04 | 0.09 | 98% | 1.983 s | 0.06 |
| blocks | node | 6 | 542.06 ms (482.19 ms–583.96 ms) | 0.23 | 0.04 | 0.12 | 97% | 1.792 s | 0.07 |
| blocks | slot | 2 | 76.30 ms (73.00 ms–86.07 ms) | 1.64 | 0.82 | 0.84 | 11% | 58.47 ms | 2.00 |
| blocks | slot | 3 | 308.58 ms (305.34 ms–387.66 ms) | 0.41 | 0.14 | 0.21 | 93% | 214.34 ms | 0.55 |
| blocks | slot | 4 | 203.83 ms (179.56 ms–217.47 ms) | 0.61 | 0.15 | 0.31 | 91% | 139.21 ms | 0.84 |
| blocks | slot | 5 | 330.85 ms (286.88 ms–331.48 ms) | 0.38 | 0.08 | 0.19 | 96% | 676.14 ms | 0.17 |
| blocks | slot | 6 | 264.31 ms (200.16 ms–301.11 ms) | 0.47 | 0.08 | 0.24 | 94% | 473.76 ms | 0.25 |
| strips | node | 2 | 360.73 ms (350.34 ms–381.70 ms) | 0.35 | 0.17 | 0.18 | 90% | 233.86 ms | 0.50 |
| strips | node | 3 | 384.13 ms (301.32 ms–440.70 ms) | 0.33 | 0.11 | 0.17 | 95% | 691.83 ms | 0.17 |
| strips | node | 4 | 470.48 ms (446.61 ms–499.48 ms) | 0.27 | 0.07 | 0.14 | 97% | 1.989 s | 0.06 |
| strips | node | 5 | 661.69 ms (611.41 ms–697.30 ms) | 0.19 | 0.04 | 0.10 | 98% | 1.984 s | 0.06 |
| strips | node | 6 | 632.65 ms (561.57 ms–780.96 ms) | 0.20 | 0.03 | 0.10 | 98% | 1.980 s | 0.06 |
| strips | slot | 2 | 74.67 ms (69.98 ms–81.82 ms) | 1.68 | 0.84 | 0.86 | 7% | 58.47 ms | 2.00 |
| strips | slot | 3 | 342.40 ms (280.99 ms–366.01 ms) | 0.37 | 0.12 | 0.19 | 92% | 214.38 ms | 0.55 |
| strips | slot | 4 | 304.78 ms (268.10 ms–348.48 ms) | 0.41 | 0.10 | 0.21 | 95% | 204.63 ms | 0.57 |
| strips | slot | 5 | 323.37 ms (323.24 ms–343.24 ms) | 0.39 | 0.08 | 0.20 | 95% | 676.24 ms | 0.17 |
| strips | slot | 6 | 317.89 ms (271.32 ms–321.90 ms) | 0.39 | 0.07 | 0.20 | 96% | 672.35 ms | 0.17 |

### Does spreading across VMs pay off?

Best measured speed-up per N (any decomposition and placement, P > 1):

| N | best S | best S_fast | at | model best S |
|---|---|---|---|---|
| 512 | 1.74 | 0.58 | blocks, map-by slot, P=2 | 1.98 (blocks, slot, P=2) |
| 2048 | 1.48 | 0.86 | strips, map-by slot, P=2 | 2.00 (blocks, slot, P=2) |
| 8192 | 1.68 | 0.86 | strips, map-by slot, P=2 | 2.00 (blocks, slot, P=2) |

## 3. Weak scaling

Work per process fixed; ideally the time per step stays constant.

| decomp | P | N | time/step (spread) | comm share | model |
|---|---|---|---|---|---|
| serial | 1 | 1025 | 1.85 ms (1.32 ms–2.00 ms) | 0% | 1.85 ms |
| strips | 2 | 1449 | 84.15 ms (76.69 ms–88.27 ms) | 98% | 69.58 ms |
| strips | 3 | 1775 | 89.49 ms (82.63 ms–101.43 ms) | 99% | 343.67 ms |
| strips | 4 | 2049 | 152.84 ms (139.67 ms–166.76 ms) | 99% | 581.77 ms |
| strips | 5 | 2291 | 187.89 ms (160.22 ms–247.25 ms) | 99% | 636.14 ms |
| strips | 6 | 2509 | 170.94 ms (156.42 ms–187.86 ms) | 99% | 685.13 ms |
| blocks | 2 | 1449 | 85.98 ms (83.01 ms–110.18 ms) | 98% | 69.54 ms |
| blocks | 3 | 1775 | 83.16 ms (82.90 ms–85.08 ms) | 99% | 343.57 ms |
| blocks | 4 | 2049 | 211.92 ms (184.97 ms–212.17 ms) | 99% | 700.59 ms |
| blocks | 5 | 2291 | 147.06 ms (119.69 ms–203.65 ms) | 99% | 635.70 ms |
| blocks | 6 | 2509 | 202.54 ms (202.20 ms–214.27 ms) | 99% | 710.03 ms |

## Prediction for a real cluster

Same compute rate, t_s = 1.0 µs, bandwidth 10 GB/s (model only).

| N | decomp | S at P=6 | S at P=64 |
|---|---|---|---|
| 512 | strips | 5.66 | 41.8 |
| 512 | blocks | 5.51 | 35.7 |
| 2048 | strips | 5.96 | 61.2 |
| 2048 | blocks | 5.96 | 61.4 |
| 8192 | strips | 5.99 | 63.5 |
| 8192 | blocks | 5.99 | 63.7 |

