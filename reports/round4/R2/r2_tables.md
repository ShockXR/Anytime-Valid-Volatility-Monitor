### T1. Oracle thresholds and realised false-alarm rates on fresh null paths (20,000 paths per law)

| Law | Method | Oracle threshold (alarm when statistic >= this) | Proved/nominal threshold | False alarms within 390 steps [95% CI] | Within 1,950 steps |
|---|---|---|---|---|---|
| Gaussian | README (plug-in, band 0.25) | 349.6 | 15,600 (= 2A) | 4.7% [4.4, 5.0] | 35.5% [34.9, 36.2] |
| Gaussian | Plug-in, band 0 | 6456 | 15,600 (= 2A) | 4.7% [4.4, 5.0] | 48.4% [47.7, 49.1] |
| Gaussian | Confidence-sequence baseline | 506.4 | 15,600 (= 2A) | 4.3% [4.1, 4.6] | 38.2% [37.6, 38.9] |
| Gaussian | Known baseline (theorem) | 4172 | 15,600 (= 2A) | 5.3% [5.0, 5.6] | 49.3% [48.6, 50.0] |
| Gaussian | Gaussian LR | 2.123e+04 | 15,600 (= 2A) | 4.9% [4.6, 5.2] | 25.3% [24.7, 25.9] |
| Gaussian | Page CUSUM (up only) | 32.91 | none | 5.4% [5.1, 5.7] | 20.9% [20.4, 21.5] |
| Gaussian | Repeated 60-step z-test | 3.773 | none | 5.2% [4.9, 5.5] | 20.8% [20.2, 21.4] |
| Student t3 | README (plug-in, band 0.25) | 471.1 | 15,600 (= 2A) | 4.9% [4.6, 5.2] | 30.0% [29.4, 30.7] |
| Student t3 | Plug-in, band 0 | 8212 | 15,600 (= 2A) | 4.9% [4.6, 5.2] | 42.4% [41.7, 43.0] |
| Student t3 | Confidence-sequence baseline | 662.1 | 15,600 (= 2A) | 4.8% [4.5, 5.1] | 29.7% [29.0, 30.3] |
| Student t3 | Known baseline (theorem) | 5781 | 15,600 (= 2A) | 4.8% [4.6, 5.1] | 38.0% [37.3, 38.7] |
| Student t3 | Gaussian LR | 2.83e+122 | 15,600 (= 2A) | 4.8% [4.5, 5.1] | 80.2% [79.6, 80.7] |
| Student t3 | Page CUSUM (up only) | 936.6 | none | 4.8% [4.5, 5.1] | 97.0% [96.7, 97.2] |
| Student t3 | Repeated 60-step z-test | 57.76 | none | 4.9% [4.6, 5.2] | 22.9% [22.3, 23.5] |
| Gaussian + jumps | README (plug-in, band 0.25) | 360.4 | 15,600 (= 2A) | 4.9% [4.6, 5.2] | 35.6% [34.9, 36.2] |
| Gaussian + jumps | Plug-in, band 0 | 6575 | 15,600 (= 2A) | 5.0% [4.7, 5.3] | 47.3% [46.6, 48.0] |
| Gaussian + jumps | Confidence-sequence baseline | 503.6 | 15,600 (= 2A) | 4.7% [4.4, 5.0] | 37.9% [37.2, 38.5] |
| Gaussian + jumps | Known baseline (theorem) | 4362 | 15,600 (= 2A) | 5.0% [4.7, 5.3] | 46.9% [46.2, 47.6] |
| Gaussian + jumps | Gaussian LR | 8.42e+58 | 15,600 (= 2A) | 5.1% [4.8, 5.4] | 26.8% [26.2, 27.5] |
| Gaussian + jumps | Page CUSUM (up only) | 412.8 | none | 5.1% [4.8, 5.4] | 38.6% [38.0, 39.3] |
| Gaussian + jumps | Repeated 60-step z-test | 33 | none | 5.0% [4.7, 5.3] | 23.4% [22.8, 24.0] |

### T2-Gaussian. Oracle tuning (5% per 390 steps on this law): median delay (IQR), detected within 1,950 steps (2,000 paths per cell)

| Method | variance x2 | variance x1.5 | variance x0.5 | variance x2/3 |
|---|---|---|---|---|
| README (plug-in, band 0.25) | 37 (24-58), 95% | 99 (54-165), 96% | 73 (58-90), 95% | 134 (97-193), 95% |
| Plug-in, band 0 | 43 (29-63), 94% | 100 (59-151), 96% | 93 (75-112), 95% | 147 (112-188), 95% |
| Confidence-sequence baseline | 44 (28-68), 95% | 122 (66-212), 96% | 70 (54-86), 96% | 120 (84-169), 95% |
| Known baseline (theorem) | 41 (27-57), 95% | 86 (53-133), 95% | 84 (68-101), 95% | 127 (97-161), 95% |
| Gaussian LR | 43 (26-65), 94% | 116 (65-191), 96% | 66 (49-87), 95% | 144 (92-225), 95% |
| Page CUSUM (up only) | 38 (24-56), 94% | 92 (53-158), 95% | n/a | n/a |
| Repeated 60-step z-test | 42 (29-56), 94% | 92 (54-182), 95% | 523 (248-976), 75% | 990 (505-1428), 7% |

### T2-Student t3. Oracle tuning (5% per 390 steps on this law): median delay (IQR), detected within 1,950 steps (2,000 paths per cell)

| Method | variance x2 | variance x1.5 | variance x0.5 | variance x2/3 |
|---|---|---|---|---|
| README (plug-in, band 0.25) | 80 (48-130), 95% | 226 (116-416), 95% | 122 (80-182), 95% | 296 (168-513), 95% |
| Plug-in, band 0 | 79 (51-119), 96% | 189 (107-301), 95% | 122 (86-170), 95% | 240 (152-366), 95% |
| Confidence-sequence baseline | 91 (55-147), 96% | 270 (133-484), 94% | 116 (75-173), 95% | 262 (151-454), 95% |
| Known baseline (theorem) | 74 (48-106), 94% | 171 (102-259), 95% | 113 (81-155), 95% | 218 (140-321), 95% |
| Gaussian LR | 249 (167-343), 95% | 422 (274-583), 95% | 936 (464-1555), 7% | 1092 (673-1568), 17% |
| Page CUSUM (up only) | 203 (134-269), 94% | 305 (200-429), 95% | n/a | n/a |
| Repeated 60-step z-test | 750 (341-1228), 67% | 856 (391-1343), 42% | 1015 (580-1569), 8% | 1026 (532-1515), 10% |

### T2-Gaussian + jumps. Oracle tuning (5% per 390 steps on this law): median delay (IQR), detected within 1,950 steps (2,000 paths per cell)

| Method | variance x2 | variance x1.5 | variance x0.5 | variance x2/3 |
|---|---|---|---|---|
| README (plug-in, band 0.25) | 39 (24-61), 95% | 102 (58-174), 95% | 77 (60-100), 95% | 150 (103-225), 94% |
| Plug-in, band 0 | 45 (30-65), 95% | 103 (63-157), 95% | 99 (77-124), 95% | 160 (116-219), 94% |
| Confidence-sequence baseline | 46 (28-72), 94% | 130 (71-224), 95% | 72 (54-94), 95% | 128 (88-191), 95% |
| Known baseline (theorem) | 42 (28-60), 95% | 92 (57-139), 94% | 89 (69-110), 94% | 140 (103-186), 94% |
| Gaussian LR | 554 (384-736), 95% | 1079 (662-1439), 78% | 1202 (668-1583), 32% | 936 (414-1384), 22% |
| Page CUSUM (up only) | 354 (234-469), 95% | 671 (415-977), 92% | n/a | n/a |
| Repeated 60-step z-test | 888 (444-1352), 35% | 886 (425-1402), 28% | 878 (373-1358), 19% | 939 (418-1397), 20% |

### T3. Where the README setting's slowness comes from (median delay in steps)

| Law | Variance | README at 2A | README at oracle | Band 0 at 2A | Band 0 at oracle | Known baseline at oracle | Best non-e competitor at oracle |
|---|---|---|---|---|---|---|---|
| Gaussian | x2 | 114 | 37 | 51 | 43 | 41 | 38 (Page CUSUM (up only)) |
| Gaussian | x1.5 | 515 | 99 | 119 | 100 | 86 | 92 (Page CUSUM (up only)) |
| Gaussian | x0.5 | 265 | 73 | 113 | 93 | 84 | 66 (Gaussian LR) |
| Gaussian | x2/3 | 860 | 134 | 179 | 147 | 127 | 144 (Gaussian LR) |
| Student t3 | x2 | 273 | 80 | 90 | 79 | 74 | 203 (Page CUSUM (up only)) |
| Student t3 | x1.5 | 966 | 226 | 221 | 189 | 171 | 305 (Page CUSUM (up only)) |
| Student t3 | x0.5 | 518 | 122 | 138 | 122 | 113 | 936 (Gaussian LR) |
| Student t3 | x2/3 | 1114 | 296 | 284 | 240 | 218 | 1026 (Repeated 60-step z-test) |
| Gaussian + jumps | x2 | 118 | 39 | 52 | 45 | 42 | 354 (Page CUSUM (up only)) |
| Gaussian + jumps | x1.5 | 559 | 102 | 123 | 103 | 92 | 671 (Page CUSUM (up only)) |
| Gaussian + jumps | x0.5 | 308 | 77 | 119 | 99 | 89 | 878 (Repeated 60-step z-test) |
| Gaussian + jumps | x2/3 | 979 | 150 | 202 | 160 | 140 | 936 (Gaussian LR) |

### T4. Practical tuning: competitors tuned on Gaussian data, e-detectors at the proved 2A

| Law | Method | Threshold used | False alarms within 390 steps [95% CI] | variance x2 | variance x1.5 | variance x0.5 | variance x2/3 |
|---|---|---|---|---|---|---|---|
| Student t3 | README (plug-in, band 0.25) | 2A (proved) | 0.0% [0.0, 0.0] | 273 (170-421), 100% | 966 (572-1367), 45% | 518 (328-776), 98% | 1114 (710-1551), 19% |
| Student t3 | Plug-in, band 0 | 2A (proved) | 2.6% [2.4, 2.9] | 90 (59-133), 98% | 221 (132-353), 97% | 138 (101-193), 97% | 284 (182-435), 97% |
| Student t3 | Confidence-sequence baseline | 2A (proved) | 0.0% [0.0, 0.0] | 254 (155-388), 100% | 891 (537-1291), 53% | 347 (229-530), 100% | 927 (573-1302), 54% |
| Student t3 | Known baseline (theorem) | 2A (proved) | 1.6% [1.5, 1.8] | 91 (61-128), 99% | 221 (140-335), 98% | 141 (104-188), 98% | 286 (189-417), 98% |
| Student t3 | Gaussian LR | tuned on Gaussian | 99.9% [99.9, 100.0] | 0.1% detected | 0.1% detected | 0.1% detected | 0.1% detected |
| Student t3 | Page CUSUM (up only) | tuned on Gaussian | 99.9% [99.9, 100.0] | 0.0% detected | 0.1% detected | n/a | n/a |
| Student t3 | Repeated 60-step z-test | tuned on Gaussian | 99.9% [99.9, 100.0] | 0.1% detected | 0.1% detected | 0.1% detected | 0.1% detected |
| Gaussian + jumps | README (plug-in, band 0.25) | 2A (proved) | 0.0% [0.0, 0.0] | 118 (77-172), 100% | 559 (336-900), 93% | 308 (240-403), 100% | 979 (652-1374), 64% |
| Gaussian + jumps | Plug-in, band 0 | 2A (proved) | 2.0% [1.9, 2.3] | 52 (36-75), 98% | 123 (77-188), 98% | 119 (96-149), 98% | 202 (145-271), 98% |
| Gaussian + jumps | Confidence-sequence baseline | 2A (proved) | 0.0% [0.0, 0.0] | 119 (76-176), 100% | 556 (323-908), 91% | 220 (172-285), 100% | 562 (375-874), 95% |
| Gaussian + jumps | Known baseline (theorem) | 2A (proved) | 1.0% [0.9, 1.2] | 53 (37-75), 99% | 125 (80-183), 99% | 119 (97-147), 99% | 198 (150-257), 99% |
| Gaussian + jumps | Gaussian LR | tuned on Gaussian | 66.6% [66.0, 67.3] | 38 (22-60), 32% | 88 (47-144), 33% | 60 (41-81), 33% | 110 (63-178), 34% |
| Gaussian + jumps | Page CUSUM (up only) | tuned on Gaussian | 65.4% [64.7, 66.1] | 34 (21-51), 33% | 75 (42-126), 35% | n/a | n/a |
| Gaussian + jumps | Repeated 60-step z-test | tuned on Gaussian | 62.6% [62.0, 63.3] | 39 (25-52), 36% | 71 (43-135), 37% | 260 (111-517), 36% | 359 (139-719), 37% |

### T5. Extra: e-detectors and GLR also tuned on Gaussian data (5% per 390 steps there), used on the other laws

| Law | Method | False alarms within 390 steps [95% CI] | variance x2 | variance x1.5 | variance x0.5 | variance x2/3 |
|---|---|---|---|---|---|---|
| Student t3 | README (plug-in, band 0.25) | 9.8% [9.4, 10.2] | 68 (38-108), 90% | 170 (86-308), 89% | 97 (63-142), 90% | 216 (123-361), 90% |
| Student t3 | Plug-in, band 0 | 6.3% [6.0, 6.6] | 75 (48-112), 94% | 178 (101-283), 94% | 117 (81-160), 93% | 224 (139-341), 94% |
| Student t3 | Confidence-sequence baseline | 8.3% [7.9, 8.7] | 79 (46-127), 92% | 213 (112-397), 91% | 98 (64-146), 91% | 214 (122-356), 92% |
| Student t3 | Known baseline (theorem) | 7.1% [6.7, 7.4] | 69 (44-99), 92% | 153 (89-232), 93% | 102 (73-141), 93% | 194 (124-289), 92% |
| Gaussian + jumps | README (plug-in, band 0.25) | 5.4% [5.1, 5.7] | 39 (23-60), 95% | 100 (56-172), 95% | 75 (58-99), 94% | 144 (98-216), 94% |
| Gaussian + jumps | Plug-in, band 0 | 5.1% [4.8, 5.5] | 45 (29-65), 95% | 103 (62-156), 95% | 98 (77-124), 95% | 159 (116-218), 94% |
| Gaussian + jumps | Confidence-sequence baseline | 4.7% [4.4, 5.0] | 46 (29-72), 94% | 130 (71-225), 95% | 72 (54-95), 95% | 129 (89-194), 95% |
| Gaussian + jumps | Known baseline (theorem) | 5.3% [5.0, 5.6] | 42 (27-59), 94% | 91 (57-137), 94% | 88 (68-109), 94% | 138 (101-185), 94% |

### T6. Nominal thresholds for reference (as in B2): e-detectors and GLR at 2A, on each law

| Law | Method | False alarms within 390 steps [95% CI] | variance x2 | variance x1.5 | variance x0.5 | variance x2/3 |
|---|---|---|---|---|---|---|
| Gaussian | README (plug-in, band 0.25) | 0.0% [0.0, 0.0] | 114 (77-164), 100% | 515 (314-799), 96% | 265 (226-311), 100% | 860 (611-1245), 78% |
| Gaussian | Plug-in, band 0 | 1.9% [1.7, 2.1] | 51 (35-72), 98% | 119 (75-183), 98% | 113 (94-132), 98% | 179 (139-227), 98% |
| Gaussian | Confidence-sequence baseline | 0.0% [0.0, 0.0] | 114 (76-164), 100% | 504 (302-794), 95% | 197 (165-235), 100% | 454 (326-631), 99% |
| Gaussian | Known baseline (theorem) | 1.0% [0.9, 1.1] | 52 (36-71), 99% | 120 (79-176), 99% | 113 (96-130), 99% | 176 (142-220), 99% |
| Gaussian | Gaussian LR | 6.2% [5.9, 6.6] | 41 (25-63), 93% | 109 (60-183), 95% | 62 (46-83), 94% | 135 (87-210), 93% |
| Student t3 | README (plug-in, band 0.25) | 0.0% [0.0, 0.0] | 273 (170-421), 100% | 966 (572-1367), 45% | 518 (328-776), 98% | 1114 (710-1551), 19% |
| Student t3 | Plug-in, band 0 | 2.6% [2.4, 2.9] | 90 (59-133), 98% | 221 (132-353), 97% | 138 (101-193), 97% | 284 (182-435), 97% |
| Student t3 | Confidence-sequence baseline | 0.0% [0.0, 0.0] | 254 (155-388), 100% | 891 (537-1291), 53% | 347 (229-530), 100% | 927 (573-1302), 54% |
| Student t3 | Known baseline (theorem) | 1.6% [1.5, 1.8] | 91 (61-128), 99% | 221 (140-335), 98% | 141 (104-188), 98% | 286 (189-417), 98% |
| Student t3 | Gaussian LR | 100.0% [99.9, 100.0] | 0.1% detected | 0.1% detected | 0.0% detected | 0.1% detected |
| Gaussian + jumps | README (plug-in, band 0.25) | 0.0% [0.0, 0.0] | 118 (77-172), 100% | 559 (336-900), 93% | 308 (240-403), 100% | 979 (652-1374), 64% |
| Gaussian + jumps | Plug-in, band 0 | 2.0% [1.9, 2.3] | 52 (36-75), 98% | 123 (77-188), 98% | 119 (96-149), 98% | 202 (145-271), 98% |
| Gaussian + jumps | Confidence-sequence baseline | 0.0% [0.0, 0.0] | 119 (76-176), 100% | 556 (323-908), 91% | 220 (172-285), 100% | 562 (375-874), 95% |
| Gaussian + jumps | Known baseline (theorem) | 1.0% [0.9, 1.2] | 53 (37-75), 99% | 125 (80-183), 99% | 119 (97-147), 99% | 198 (150-257), 99% |
| Gaussian + jumps | Gaussian LR | 67.5% [66.9, 68.2] | 36 (21-58), 31% | 84 (42-138), 33% | 58 (38-78), 32% | 105 (60-167), 33% |

### T7. Matched over five days instead of one: threshold that gives 25% false alarms within 1,950 steps on this law

| Law | Method | False alarms within 390 / 1,950 steps | variance x2 | variance x1.5 | variance x0.5 | variance x2/3 |
|---|---|---|---|---|---|---|
| Gaussian | README (plug-in, band 0.25) | 3.4% / 25.2% | 40 (25-61) | 108 (58-180) | 80 (65-99) | 153 (110-219) |
| Gaussian | Plug-in, band 0 | 1.6% / 24.9% | 52 (36-73) | 123 (77-188) | 116 (97-135) | 183 (144-233) |
| Gaussian | Confidence-sequence baseline | 2.6% / 25.0% | 49 (30-74) | 140 (77-245) | 78 (61-96) | 136 (97-195) |
| Gaussian | Known baseline (theorem) | 2.2% / 25.8% | 46 (32-64) | 105 (66-154) | 99 (82-116) | 151 (121-191) |
| Gaussian | Gaussian LR | 5.0% / 25.4% | 43 (26-65) | 115 (65-191) | 65 (49-86) | 144 (92-224) |
| Gaussian | Page CUSUM (up only) | 6.7% / 25.3% | 36 (23-54) | 85 (48-146) | n/a | n/a |
| Gaussian | Repeated 60-step z-test | 6.4% / 25.0% | 41 (28-55) | 84 (52-164) | 356 (160-694) | 885 (405-1358) |
| Student t3 | README (plug-in, band 0.25) | 4.0% / 25.2% | 85 (50-136) | 240 (123-446) | 130 (85-192) | 323 (179-554) |
| Student t3 | Plug-in, band 0 | 2.2% / 25.1% | 93 (61-137) | 230 (138-368) | 144 (105-201) | 294 (191-453) |
| Student t3 | Confidence-sequence baseline | 3.9% / 25.1% | 96 (58-153) | 288 (147-524) | 123 (79-185) | 278 (160-489) |
| Student t3 | Known baseline (theorem) | 3.0% / 25.2% | 81 (53-118) | 193 (118-296) | 125 (92-169) | 249 (163-364) |
| Student t3 | Gaussian LR | 1.2% / 25.0% | 535 (398-689) | 929 (653-1175) | 2.5% detected | 4.7% detected |
| Student t3 | Page CUSUM (up only) | 0.6% / 24.9% | 662 (530-805) | 1011 (788-1214) | n/a | n/a |
| Student t3 | Repeated 60-step z-test | 5.4% / 24.9% | 695 (322-1196) | 840 (396-1332) | 1000 (624-1546) | 1015 (504-1511) |
| Gaussian + jumps | README (plug-in, band 0.25) | 3.4% / 24.7% | 42 (26-66) | 112 (63-188) | 85 (66-111) | 176 (118-259) |
| Gaussian + jumps | Plug-in, band 0 | 1.8% / 25.6% | 53 (36-76) | 126 (78-191) | 121 (97-151) | 206 (150-276) |
| Gaussian + jumps | Confidence-sequence baseline | 2.9% / 24.8% | 50 (31-80) | 147 (84-262) | 81 (61-106) | 150 (104-225) |
| Gaussian + jumps | Known baseline (theorem) | 2.3% / 25.2% | 47 (32-68) | 108 (68-162) | 104 (83-129) | 170 (129-223) |
| Gaussian + jumps | Gaussian LR | 4.8% / 25.5% | 570 (397-754) | 1102 (679-1457) | 1199 (654-1564) | 932 (414-1390) |
| Gaussian + jumps | Page CUSUM (up only) | 2.3% / 25.2% | 470 (330-605) | 841 (560-1219) | n/a | n/a |
| Gaussian + jumps | Repeated 60-step z-test | 5.5% / 25.5% | 888 (440-1366) | 888 (422-1392) | 860 (369-1344) | 940 (423-1378) |
