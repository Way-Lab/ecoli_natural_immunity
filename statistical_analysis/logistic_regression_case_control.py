"""
Matched Case-Control Analysis of Maternal Antibodies and Neonatal Sepsis
==========================================================================

Unified script for analyzing multiple antibody assays:
- THP-1 (opsonization activity)
- Mix8 (IgG EcEPT)
- HL60 (opsonization activity)
- OmpA (IgG EcEPT)

Key Features:
- Handles 1:3 matched case-control design
- Conditional logistic regression (accounts for matching)
- Calculates protective thresholds
- Bootstrap confidence intervals


Usage:
    python matched_case_control_antibody_analysis.py                  # Run all assays
    python matched_case_control_antibody_analysis.py --assay THP1     # Run single assay
    python matched_case_control_antibody_analysis.py --assay Mix8 HL60  # Run multiple assays
    python matched_case_control_antibody_analysis.py --list           # List available assays
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from io import StringIO
import statsmodels.api as sm
from statsmodels.discrete.conditional_models import ConditionalLogit
import warnings
import argparse
from dataclasses import dataclass
from typing import Optional

warnings.filterwarnings('ignore')

# =============================================================================
# CONFIGURATION CLASSES
# =============================================================================

@dataclass
class AssayConfig:
    """Configuration for each antibody assay analysis."""
    name: str
    data_str: str
    output_figure: str
    output_summary: str
    y_axis_label: str
    y_axis_bottom: float  # 10^x for log scale bottom
    x_axis_max: float     # Upper limit for probability curve panel
    use_standard_logistic_for_viz: bool = False  # HL60 uses different method

# =============================================================================
# ASSAY DATA AND CONFIGURATIONS
# =============================================================================

MIX8_DATA = """Case-Mix8	control 1	control 2	control 3
14874.8639	5974.17582	101256.987	9590.6049
17166.7467	14630.3933	15923.1138	13737.6675
21757.801	10422.4303	27606.0964	2199.9307
313.833283	77788.0707	6269.64571	5170.74356
2127.5364	313.833283	10386.4707	902.632569
1742.98641	1669.66368	21116.8921	9359.48602
9220.63859	8958.37009	13737.6675	8703.08049
3896.90663	7418.00292	6414.7096	8697.11074
5112.44913	3201.30461	46531.1729	30855.8134
3148.61094	11806.74	3374.18585	1106.83275
1424.2263	6972.47555	6485.76188	16595.143
13849.8797	12012.2953	4680.36799	4984.81448
1384.98797	1424.2263	3451.59634	3048.98298
4979.02148	8879.95981	7500.92806	902.632569
1168.89634	750.092806	887.995981	313.833283
3438.47412	79560.9463	2666.44807	10245.9209
1084.63411	49511.6934	2549.79949	5457.97799
10595.0241	9264.36975	9316.14137	29534.8126
16137.0926	17067.4618	17067.4618	24131.3984
2707.88016	8440.45016	441.967668	4302.48683
935.807432	8440.45016	1089.56135	2707.88016
11758.9393	24952.8485	441.967668	49511.6934
327.757479	356.536601	20789.5084	14932.1169
309.029543	4951.16934	777.887244	481.59571
313.833283	5982.66934		
309.029543	481.59571		
1113.933	7301.651	1578.673	1799.221
768.7318	902.1879	2722.82	530.0315
755.5082	6655.959	3790.976	2299.412
398.9596	5515.014	2645.916	1234.216
550.8933	1464.721	1479.898	733.4319
235.7199	7518.695	2670.43	942.952
751.5659	3204.696	17837.4	4440.956
978.1437	5313.144	1485.36	2029.521
2054.765	4377.684	1557.701	915.4313
627.4738	321.4153	583.7767	1430.29
330.2284	4799.962	1604.884	3640.231
1526.117	7913.169	29240.1	9800.078
4049.324	7545.375	4185.052	2679.341
125.3556	629.1067	186.1947	1439.802
125.3556	3382.276	794.3315	394.0294
306.6076	2384.25	958.1306	742.2476
125.3556	2441.678	1160.534	6526.849
125.3556	390.4786	1481.016	200.3712
125.3556	598.3896	3078.715	125.3556
125.3556	1270.334	8077.957	933.7587
414.0601	2193.561	3823.669	15908.14
797.1428	208.3054	1036.631	3187.916
780.898352	1039.53596	899.80342	759.231639
125.3556	198.9137	1101.885	3129.441
125.3556	818.248284	613.609022	836.432969
1212.99708	14548.5277	3353.27992	4964.48167
269.544792	1043.03506	220.51217	1865.069022
2507.112	67645.52	15886.63	7880.588
125.3555672	1035.213258	523.964016	125.3555672
1228.584378	25098.0384	10620.24718	6836.51624
4588.48212	1747.207868	1460.254756	1671.588842
2642.28028	2726.23284	4245.61226	928.683848
14945.09132	6848.37382	5733.17914	3418.64244
3243.18818	5437.03576	7843.56558	8694.63498
1063.33141	1027.995024	1478.52886	1757.103728
607.27443	2934.97642	2228.99548	1357.317834
2874.17074	2563.71172	3828.52318	6802.02776
1493.075546	1842.76015	1792.30234	2336.30488
1264.173148	3789.3938	4571.59088	2520.79252
2161.483	7336.87246	6332.39817	4847.09647
2025.161	20682.2684	7804.52324	13809.491
2096.558	4021.12125	6842.4877	5974.7188
1636.206	26960.8153	76728.8848	68179.097
1739.314	4402.72239	2497.62816	1746.91799
841.7201	1157.01151	9080.86736	6424.69566
1207.322	16113.6577	26539.3092	14263.6961
18740.43	16934.0019	45949.9015	18497.8487
1267.089	7943.09436	24033.8561	421.180514
571.1765	7288.75225	10609.1473	12610.7712
394.4413	7688.17772	9325.07671	7922.89275
502.0907	12303.0483	674.45231	18648.8411
198.7618	1305.88354	1221.68959	1015.62831
229.1332	759.068253	349.966106	125.355567
196.0579	2605.4231	433.691204	4423.58864
130.6058	735.225087	587.658826	304.415499
125.3556	6662.55169	207.460459	5777.89684
125.3556	174.605896	377.944664	573.870062
125.3556	439.657826	2292.1512	243.76872
210.4719	3652.13698	3888.49037	28425.4441
266.2733	1610.37883	360.306831	626.362548
500.4355	963.498196	1484.17176	87307.8245
4585.238	1312.47712	13130.1251	82549.8611
343.4899	29784.9155	16929.1156	12529.5349
291.0112	5388.22946	39150.9405	5982.38372
191.7986	2093.3663	3487.21246	1718.73723
145.705927	2926.5406	4417.99207	14759.3421
7240.50388	42108.2884	46240.8124	8953.5581
319.200764	9538.31395	1592.49677	948.166788
303.500271	2392.52469	1624.80624	842.671982
119.393807	1126.41123	401.597955	508.800775
185.403555	3384.54111	2290.68343	600.124673
333.336023	5385.41675	623.65471	484.500661
125.3556	3514.53192	617.968534	389.361874
817.952573	14400.3228	1817.64786	1012.19252
"""

OMPA_DATA = """Case-OmpA Loop	control 1	control 2	control 3
696.6426575	3010.357059	1202.858063	1840.825936
145.3964974	1039.919833	2062.083732	1744.810087
204.7481229	1408.463733	3311.261423	174.7671906
77.05293468	652.6935716	327.2885588	77.05293468
664.6965502	4323.686784	1725.925768	6619.262181
883.0736768	220.141416	327.1961358	315.5517578
527.7173978	136.0619735	301.4133193	2024.994396
234.8762131	4280.48632	2981.104638	4778.529468
1227.602395	5049.010982	838.781499	2563.978057
285.781357	121.0349442	944.4369787	1037.056801
312.1892578	1382.864514	2526.118321	206.291253
283.3314403	1661.656659	579.676496	224.6967629
499.6192698	690.8147276	5390.027058	3199.880471
207.8936899	257.6426632	10326.95497	1348.868974
82.81162915	188.6007457	1165.395476	311.4413493
77.05293468	206.7720794	2089.850897	1699.87577
66.5648591	791.2640985	1778.258939	226.2097354
319.6259692	606.1490832	179.330591	948.937191
385.862674	726.365842	861.3901881	69.23885394
375.0022662	627.0814792	822.1316531	692.2916274
102.8759811	152.5262295	860.169025	2080.127724
155.2791874	104.7644684	294.82138	1277.591762
154.7382136	77.05293468	2607.942236	1421.59185
77.05293468	1003.327526	603.1346709	51.44506805
136.0196272	679.8694013		
155.4508214	1230.188441		
319.9942	694.5039	1323.479	7098.704
80.42745	340.6539	108.1904	9847.488
75.83385	303.7703	188.5462	155.7038
261.7506	1107.43	212.4758	1869.193
205.886	233.4309	1978.114	223.2642
43.97227	204.2698	23.16707	37.64189
95.37095	88.38398	156.9419	374.8215
385.862	637.3182	30.30328	303.9923
451.1955	353.1754	5120.823	1438.471
626.175	560.8157	3045.007	638.3704
287.8914	1017.068	650.2168	411.4742
415.5476	1957.464	3694.005	2361.389
956.6345	641.6826	2959.402	469.6529
673.2549	434.8593	2886.649	1442.134
1178.846	5755.757	5569.955	853.9744
879.1198	588.878	3428.353	817.9824
86.5773	426.2104	75.69619	217.4964
457.5104	2786.192	840.8993	141.8061
105.9354	133.3825	30.19493	506.7911
95.46123	1646.3	177.0549	115.575
246.6448	405.7934	413.7395	976.6142
661.6952	500.0276	323.3375	93.16247
44.04139	627.8097	271.2083	501.1607
48.91038	145.6696	417.364	471.673
22.73804	94.95919	10.16675	720.3375
733.1942	662.7077	315.6334	721.2653
89.98193	51.82626	1433.074	381.9545
184.8937	470.8437	367.1895	446.0219
274.9682	1020.422	875.7287	1344.538
503.8322	1811.395	3601.167	715.2805
336.7617	55.43771	115.6238	252.0048
42.59329	57.78311	55.44943	994.5265
155.9911	80.36074	652.9431	508.2862
250.4115	295.937	414.9623	2182.06
223.3958	51.10993	406.7174	1846.318
516.1735	605.8124	714.7325	272.7337
642.1894	1052.977	430.7596	440.12
445.9707	3775.437	1897.184	803.6373
66.23731	186.5886	259.0459	365.9569
94.66892	58.39404	164.3323	309.879
861.5471	1226.099	2152.388	749.0991
635.7265	171.2033	6274.678	994.113
191.557226	26.2338511	213.659235	971.40244
753.612633	944.214923	2306.69392	42.0506051
742.21022	7219.51412	8312.88208	46.8675024
2117.26452	2058.31572	1965.68462	4474.40498
234.166028	317.29826	118.075409	165.112065
548.875834	1705.8456	251.222018	1194.16613
492.439283	348.901543	1622.56786	342.135335
339.100625	225.930779	1267.05632	1231.50356
28.8043202	162.578327	107.683921	254.225418
157.602151	1081.14243	582.796842	88.4493488
471.083063	342.804639	5424.38892	3098.44734
393.212042	1625.02217	2861.8371	1699.26045
255.7769	201.276332	628.17182	18.6412105
200.470304	196.765988	1569.19812	1014.9321
772.787307	1604.31831	1781.27975	306.241815
325.634229	349.184205	299.878466	747.464808
248.846929	2104.51894	481.233691	1010.27869
94.684777	936.844208	1800.85977	4893.74276
356.738552	739.577454	2549.8103	1458.41162
1672.35272	489.231889	898.67189	785.269394
359.37937	547.455893	453.769695	7735.99885
1219.44452	1468.42736	779.073903	2569.98316
289.434996	840.678217	429.118696	390.3171
632.236328	2498.68682	689.639873	143.998429
351.423938	2131.33666	675.988826	962.842394
477.468931	480.543635	2714.42266	369.452559
216.513512	2498.4301	2085.4404	82.6271243
234.363223	568.01032	461.219698	2137.77287
378.305514	585.194895	273.985225	1723.84211
284.715854	1480.97427	2305.87239	993.126629
386.215349	1798.80167	2026.07403	1096.6775
301.510631	265.910786	378.735847	346.341398
"""

HL60_DATA = """Case-HL60	control 1	control 2	control 3
114	161	249	182
100	170	187	160
125	187	186	166
38	292	137	154
102	157	238	190
122	165	230	232
53	211	226	202
87	151	141	162
69	154	311	379
84	201	160	154
64	165	139	164
45	202	202	133
20	168	234	197
124	140	138	177
121	122	148	159
38	160	145	146
29	255	176	194
48	307	227	331
45	282	273	304
5	184	121	180
21	132	180	186
12	148	87	216
60	77	205	164
42	128	76	231
2	183		
63	81		
54	476	180	159
10	50	377	91
10	388	92	131
40	67	179	118
10	105	86	117
70	96	166	105
102	157	126	173
61	113	103	95
67	117	121	77
64	103	99	161
67	140	109	108
104	124	113	119
91	134	126	114
82	131	88	92
62	71	108	66
70	94	113	93
101	160	145	128
80	163	190	89
108	215	145	122
61	126	140	206
63	304	142	201
73	281	254	382
33	176	152	180
160	236	372	276
121	174	142	375
40	193	358	212
96	281	201	336
9	189	315	281
56	148	60	58
3	134	377	127
157	227	188	181
104	183	183	166
141	185	230	151
57	112	72	74
51	96	114	241
57	118	185	88
72	126	86	136
160	215	183	205
17	92	347	141
44	128	169	155
53	376	99	371
60	212	195	215
52	249	444	383
17	156	187	211
51	71	219	175
67	190	190	229
42	149	178	159
37	143	180	99
66	75	136	230
46	113	130	131
13	162	135	213
20	155	167	77
33	151	91	82
26	330	243	269
5	163	243	178
89	237	145	184
3	157	257	147
8	280	397	188
22	307	371	313
34	185	33	168
79	173	183	379
113	191	209	345
112	276	228	225
54	206	276	223
53	200	275	200
37	200	207	366
81	313	264	200
32	196	168	183
73	196	162	180
17	193	160	180
41	190	155	185
102	185	171	184
23	168	166	184
129	235	176	171
"""

THP1_DATA = """Case-THP1	control 1	control 2	control 3
22.33	8.96598	658.332	109.3879
60.75	504.3016	169.1158	169.1158
47.66	148.5734	221.8804	16.248
4.036	504.3016	31.41025	26.71691
10.71	1	103.6768	21.62227
10.12	14.97268	353.4242	111.3602
66.87	121.7811	328.5777	115.4137
13.72	69.44533	38.19623	39.09078
16.25	18.87882	391.2982	381.2601
15.61	98.25639	20.92546	88.2033
5.096	80.54082	51.79198	267.8273
76.22	96.5104	70.76081	60.749
11.89	14.97268	19.55385	14.34469
51.79	80.54082	94.79354	10.12041
12.5	14.97268	17.54977	3.513225
35.58	412.6595	16.248	61.93916
11.3	86.62171	217.4571	30.60566
50.74	101.8388	121.7811	801.1796
69.45	74.82606	362.3451	371.6167
14.97	193.2144	13.72302	17.54977
13.72	72.0958	13.72302	33.88053
94.79	115.4137	4.563172	185.9137
2.482	55.03805	353.4242	133.235
1.47	2.482466	16.248	18.21083
1.974118	25.96496		
1.470374	16.248		
5.632849	68.14893	52.85966	83.53334
6.175261	16.89553	13.10756	43.7314
20.92546	101.8388	100.0323	33.04752
16.248	461.3617	143.2495	33.04752
10.12041	119.6189	38.19623	96.5104
20.92546	634.9777	4059.248	101.8388
30.60566	70.76081	172.3207	117.4968
41.84075	208.9557	115.4137	17.54977
43.7314	93.10509	48.67167	83.53334
47.65826	94.79354	105.547	32.22415
27.47727	47.65826	267.8273	64.37017
44.69458	36.4391	79.07913	107.4504
79.07913	21.62227	119.6189	86.62171
50.73827	98.25639	49.69826	43.7314
2.99548	109.3879	21.62227	37.3124
44.69458	66.8712	49.69826	65.61174
20.92546	401.7558	1535.998	44.69458
32.22415	93.10509	36.4391	77.63977
80.54082	73.45073	138.141	74.82606
53.94159	245.9239	25.22128	33.04752
1.470374	50.73827	10.70612	23.75819
20.23602	55.03805	59.57527	58.41764
50.73827	61.93916	148.5734	169.1158
18.21083	98.25639	89.81063	60.749
23.03851	178.957	156.994	38.19623
68.14893	82.02538	64.37017	64.37017
135.6634	267.8273	94.79354	279.8243
126.2303	1535.998	1009.349	801.1796
138.141	231.0963	148.5734	138.141
103.6768	658.332	353.4242	313.4446
64.37017	145.8837	148.5734	151.3205
103.6768	182.3939	172.3207	74.82606
448.3535	226.4251	204.8688	197.0019
64.37017	148.5734	162.92	344.8345
13.10756	44.69458	50.73827	52.85966
7.275592	133.235	105.547	123.9845
115.4137	126.2303	185.9137	240.8399
55.03805	68.14893	74.82606	89.81063
38.86764	80.46919	123.9218	74.74186
57.13634	117.4425	73.3631	66.76507
57.13634	469.4716	81.95658	115.3613
14.59685	126.1639	284.7431	31.14365
31.14365	80.46919	130.779	94.7417
68.04667	1364.496	278.5436	126.1639
8.544701	44.50012	101.7906	51.63015
35.33371	121.7216	483.3407	94.7417
77.56209	49.52723	782.0184	156.8359
31.96239	85.0017	182.1027	15.23588
11.49619	77.56209	94.7417	185.6
11.49619	85.0017	91.38934	52.70239
113.368536	634.977681	197.00192	572.728142
197.00192	1951.8779	154.12668	1259.76458
140.669137	1654.63264	226.425122	126.230322
53.9415869	572.728142	185.913659	306.263661
7.83365904	61.9391587	21.6222696	18.8788229
29.810224	93.105086	30.6056637	1062.95286
5.6328493	96.5103995	30.6056637	64.3701678
13.7230165	29.023779	42.7801978	40.9128346
21.6222696	30.6056637	178.956992	2649.12457
46.6577724	235.899358	29.023779	148.573437
29.023779	74.8260586	208.95568	683.199399
29.023779	200.885463	53.9415869	2141.52459
37.3123956	391.298232	2999.6908	213.150346
22.3265814	35.5761391	1009.34869	44.6945777
4.03576649	83.5333365	6.72281728	23.0385128
48.6716733	65.6117433	33.880528	213.150346
100.032265	836.375677	424.038733	182.393852
7.83365904	489.287826	256.548047	154.12668
1	235.899358	169.115778	10.1204092
11.8949434	88.203301	14.9726845	4.5631722
3.51322462	683.199399	154.12668	14.3446927
8.96597999	1791.70589	8.39709586	19.5538477
5.09550997	172.320657	30.6056637	2.48246568
37.3123956	1535.99754	336.557993	217.457066
"""

# Define all assay configurations
ASSAY_CONFIGS = {
    'THP1': AssayConfig(
        name='THP-1',
        data_str=THP1_DATA,
        output_figure='THP-1_analysis_figure.pdf',
        output_summary='THP-1_analysis_summary.txt',
        y_axis_label='anti-E. coli opsonization activity (THP-1 cells)',
        y_axis_bottom=-1,  # 10^-1
        x_axis_max=500,
    ),
    'Mix8': AssayConfig(
        name='Mix8',
        data_str=MIX8_DATA,
        output_figure='Mix8_antibody_analysis_figure.pdf',
        output_summary='Mix8_analysis_summary.txt',
        y_axis_label='anti-E. coli IgG EcEPT',
        y_axis_bottom=2,  # 10^2
        x_axis_max=20000,
    ),
    'HL60': AssayConfig(
        name='HL60',
        data_str=HL60_DATA,
        output_figure='HL60_opsonization_analysis_figure.pdf',
        output_summary='HL60_analysis_summary.txt',
        y_axis_label='anti-E. coli opsonization activity (HL60 cells)',
        y_axis_bottom=0,  # 10^0
        x_axis_max=500,
        use_standard_logistic_for_viz=True,
    ),
    'OmpA': AssayConfig(
        name='OmpA',
        data_str=OMPA_DATA,
        output_figure='OmpA_antibody_analysis_figure.pdf',
        output_summary='OmpA_analysis_summary.txt',
        y_axis_label='anti-OmpA IgG EcEPT',
        y_axis_bottom=1,  # 10^1
        x_axis_max=4000,
    ),
}

# Plot style constants
CASE_COLOR = '#FF2600'
CONTROL_COLOR = '#0433FF'
FIGURE_DPI = 300


# =============================================================================
# ANALYSIS FUNCTIONS
# =============================================================================

def setup_plot_style():
    """Configure matplotlib for publication-quality figures."""
    plt.style.use('seaborn-v0_8-whitegrid')
    sns.set_palette("colorblind")
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = 17


def load_and_preprocess_data(data_str: str) -> pd.DataFrame:
    """Load data from string and convert to long format."""
    df_wide = pd.read_csv(StringIO(data_str), sep='\t')
    
    records = []
    for idx, row in df_wide.iterrows():
        # Add case (first column)
        if pd.notna(row.iloc[0]):
            records.append({
                'match_id': idx + 1,
                'outcome': 1,
                'antibody': row.iloc[0]
            })
        
        # Add controls (columns 1-3)
        for i in range(1, 4):
            if i < len(row) and pd.notna(row.iloc[i]):
                records.append({
                    'match_id': idx + 1,
                    'outcome': 0,
                    'antibody': row.iloc[i]
                })
    
    df_long = pd.DataFrame(records)
    df_long = df_long.dropna()
    df_long['log_antibody'] = np.log10(df_long['antibody'])
    
    return df_long


def compute_descriptive_stats(df_long: pd.DataFrame) -> dict:
    """Compute descriptive statistics for cases and controls."""
    cases = df_long[df_long['outcome'] == 1]['antibody']
    controls = df_long[df_long['outcome'] == 0]['antibody']
    
    return {
        'cases': cases,
        'controls': controls,
        'cases_gmean': stats.gmean(cases),
        'controls_gmean': stats.gmean(controls),
        'gmr': stats.gmean(controls) / stats.gmean(cases),
    }


def compute_matched_differences(df_long: pd.DataFrame) -> np.ndarray:
    """Compute within-matched-set differences (case - mean of controls)."""
    matched_diffs = []
    for match_id in df_long['match_id'].unique():
        matched_set = df_long[df_long['match_id'] == match_id]
        case_val = matched_set[matched_set['outcome'] == 1]['log_antibody'].values
        ctrl_vals = matched_set[matched_set['outcome'] == 0]['log_antibody'].values
        
        if len(case_val) == 1 and len(ctrl_vals) > 0:
            diff = case_val[0] - np.mean(ctrl_vals)
            matched_diffs.append(diff)
    
    return np.array(matched_diffs)


def fit_conditional_logistic(df_long: pd.DataFrame):
    """Fit conditional logistic regression model."""
    model = ConditionalLogit(
        endog=df_long['outcome'],
        exog=df_long[['log_antibody']],
        groups=df_long['match_id']
    )
    return model.fit(method='bfgs', maxiter=100, disp=False)


def generate_predictions_conditional(df_long: pd.DataFrame, result, n_bootstrap: int = 1000):
    """Generate predictions using conditional logistic regression with bootstrap CI."""
    ab_range = np.logspace(
        np.log10(df_long['antibody'].min()),
        np.log10(df_long['antibody'].max()),
        300
    )
    log_ab_range = np.log10(ab_range)
    
    # Point estimates
    baseline_log_odds = result.params['log_antibody'] * np.median(df_long['log_antibody'])
    log_odds = result.params['log_antibody'] * log_ab_range
    baseline_prevalence = df_long['outcome'].mean()
    prob_pred = 1 / (1 + np.exp(-(log_odds - baseline_log_odds + 
                                  np.log(baseline_prevalence / (1 - baseline_prevalence)))))
    
    # Bootstrap
    prob_boot = np.zeros((len(ab_range), n_bootstrap))
    np.random.seed(42)
    successful_boots = 0
    
    for i in range(n_bootstrap):
        match_ids = df_long['match_id'].unique()
        boot_match_ids = np.random.choice(match_ids, size=len(match_ids), replace=True)
        
        boot_data = []
        for new_id, orig_match_id in enumerate(boot_match_ids):
            matched_set = df_long[df_long['match_id'] == orig_match_id].copy()
            matched_set['match_id'] = new_id
            boot_data.append(matched_set)
        
        df_boot = pd.concat(boot_data, ignore_index=True)
        
        try:
            model_boot = ConditionalLogit(
                endog=df_boot['outcome'],
                exog=df_boot[['log_antibody']],
                groups=df_boot['match_id']
            )
            result_boot = model_boot.fit(method='bfgs', maxiter=100, disp=False)
            
            baseline_lo_boot = result_boot.params['log_antibody'] * np.median(df_boot['log_antibody'])
            lo_boot = result_boot.params['log_antibody'] * log_ab_range
            prob_boot[:, i] = 1 / (1 + np.exp(-(lo_boot - baseline_lo_boot +
                                               np.log(df_boot['outcome'].mean() / (1 - df_boot['outcome'].mean())))))
            successful_boots += 1
        except Exception:
            prob_boot[:, i] = prob_pred
    
    prob_lower = np.percentile(prob_boot, 2.5, axis=1)
    prob_upper = np.percentile(prob_boot, 97.5, axis=1)
    
    return ab_range, prob_pred, prob_lower, prob_upper, successful_boots


def generate_predictions_standard_logistic(df_long: pd.DataFrame, n_bootstrap: int = 1000):
    """Generate predictions using standard logistic regression with bootstrap CI (for HL60)."""
    ab_range = np.logspace(
        np.log10(df_long['antibody'].min()),
        np.log10(df_long['antibody'].max()),
        300
    )
    log_ab_range = np.log10(ab_range)
    
    # Fit standard logistic regression
    X_vis = sm.add_constant(df_long['log_antibody'])
    model_vis = sm.Logit(df_long['outcome'], X_vis)
    result_vis = model_vis.fit(disp=False)
    
    X_pred = sm.add_constant(log_ab_range)
    prob_pred = result_vis.predict(X_pred)
    
    # Bootstrap
    prob_boot = np.zeros((len(ab_range), n_bootstrap))
    np.random.seed(42)
    successful_boots = 0
    
    for i in range(n_bootstrap):
        boot_idx = np.random.choice(len(df_long), size=len(df_long), replace=True)
        df_boot = df_long.iloc[boot_idx]
        
        try:
            X_boot = sm.add_constant(df_boot['log_antibody'])
            model_boot = sm.Logit(df_boot['outcome'], X_boot)
            result_boot = model_boot.fit(disp=False)
            prob_boot[:, i] = result_boot.predict(X_pred)
            successful_boots += 1
        except Exception:
            prob_boot[:, i] = prob_pred
    
    prob_lower = np.percentile(prob_boot, 2.5, axis=1)
    prob_upper = np.percentile(prob_boot, 97.5, axis=1)
    
    return ab_range, prob_pred, prob_lower, prob_upper, successful_boots


def create_figure_or_rrr(config: AssayConfig, df_long: pd.DataFrame,
                          ab_range: np.ndarray, result, n_bootstrap: int = 1000):
    """Create figure with Odds Ratio and Relative Risk Reduction curves.

    These metrics are more interpretable for case-control studies than
    absolute probability, which depends on the sampling ratio.
    """
    coef = result.params['log_antibody']
    log_ab_range = np.log10(ab_range)

    # Use minimum value as reference (low antibody = high risk baseline)
    reference_ab = df_long['antibody'].min()
    log_reference = np.log10(reference_ab)

    # Calculate Odds Ratio relative to reference
    # OR = exp(coef * (log_ab - log_reference))
    # At reference: OR = 1.0
    # As antibody increases (and coef is negative): OR decreases
    or_pred = np.exp(coef * (log_ab_range - log_reference))

    # Calculate Relative Risk Reduction
    # RRR = (1 - OR) * 100 when OR < 1 (protective)
    # Capped at 0-100% range
    rrr_pred = np.clip((1 - or_pred) * 100, 0, 100)

    # Bootstrap for confidence intervals
    np.random.seed(42)
    or_boot = np.zeros((len(ab_range), n_bootstrap))

    for i in range(n_bootstrap):
        match_ids = df_long['match_id'].unique()
        boot_match_ids = np.random.choice(match_ids, size=len(match_ids), replace=True)

        boot_data = []
        for new_id, orig_match_id in enumerate(boot_match_ids):
            matched_set = df_long[df_long['match_id'] == orig_match_id].copy()
            matched_set['match_id'] = new_id
            boot_data.append(matched_set)

        df_boot = pd.concat(boot_data, ignore_index=True)

        try:
            model_boot = ConditionalLogit(
                endog=df_boot['outcome'],
                exog=df_boot[['log_antibody']],
                groups=df_boot['match_id']
            )
            result_boot = model_boot.fit(method='bfgs', maxiter=100, disp=False)
            coef_boot = result_boot.params['log_antibody']
            or_boot[:, i] = np.exp(coef_boot * (log_ab_range - log_reference))
        except Exception:
            or_boot[:, i] = or_pred

    or_lower = np.percentile(or_boot, 2.5, axis=1)
    or_upper = np.percentile(or_boot, 97.5, axis=1)
    rrr_lower = np.clip((1 - or_upper) * 100, 0, 100)  # Note: inverted for RRR
    rrr_upper = np.clip((1 - or_lower) * 100, 0, 100)

    # Create 3-panel figure
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 8),
                                         gridspec_kw={'width_ratios': [0.5, 1, 1]})

    # Panel A: Distribution (same as original)
    np.random.seed(42)
    cases_data = df_long[df_long['outcome'] == 1]['antibody']
    controls_data = df_long[df_long['outcome'] == 0]['antibody']

    jitter_amount = 0.1
    x_controls = np.random.normal(0, jitter_amount, size=len(controls_data))
    x_cases = np.random.normal(1, jitter_amount, size=len(cases_data))

    ax1.scatter(x_controls, controls_data, alpha=0.6, s=40,
               color=CONTROL_COLOR, edgecolors='black', linewidth=0.5,
               label=f'Controls (n={len(controls_data)})')
    ax1.scatter(x_cases, cases_data, alpha=0.6, s=40,
               color=CASE_COLOR, edgecolors='black', linewidth=0.5,
               label=f'Cases (n={len(cases_data)})')

    ax1.hlines(stats.gmean(controls_data), -0.2, 0.2,
              colors='black', linewidth=3, label='Geometric mean')
    ax1.hlines(stats.gmean(cases_data), 0.8, 1.2,
              colors='black', linewidth=3)

    ax1.set_yscale('log')
    ax1.set_ylim(bottom=10**config.y_axis_bottom)
    ax1.set_xlim(-0.5, 1.5)
    ax1.set_xticks([0, 1])
    ax1.set_xticklabels(['Control', r'$\it{E. coli}$' + '\nsepsis'], fontsize=16, fontweight='bold')
    ax1.set_ylabel(config.y_axis_label, fontsize=16, fontweight='bold')
    ax1.tick_params(axis='y', labelsize=14)
    ax1.grid(True, alpha=0.3, axis='y')
    ax1.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), frameon=True, fontsize=11, ncol=1)
    ax1.set_title('A. Antibody Distribution', fontsize=18, fontweight='bold', pad=10)

    # Panel B: Odds Ratio
    ax2.fill_between(ab_range, or_lower, or_upper,
                    alpha=0.3, color='steelblue', label='95% CI', edgecolor='none')
    ax2.plot(ab_range, or_pred, 'b-', linewidth=2.5, label='Odds Ratio')
    ax2.axhline(y=1.0, color='gray', linestyle='--', linewidth=1.5, label='Reference (OR=1)')
    ax2.axvline(x=reference_ab, color='red', linestyle=':', linewidth=1.5,
                label=f'Reference: {reference_ab:.0f}')

    ax2.set_xlabel(config.y_axis_label, fontsize=16, fontweight='bold')
    ax2.set_ylabel('Odds Ratio\n(relative to minimum)', fontsize=16, fontweight='bold')
    ax2.tick_params(axis='both', labelsize=14)
    ax2.set_xlim(0, config.x_axis_max)
    ax2.set_ylim(0, 1.5)
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc='upper right', frameon=True, fontsize=11)
    ax2.set_title('B. Odds Ratio vs Antibody Level', fontsize=18, fontweight='bold', pad=10)

    # Panel C: Relative Risk Reduction
    ax3.fill_between(ab_range, rrr_lower, rrr_upper,
                    alpha=0.3, color='darkgreen', label='95% CI', edgecolor='none')
    ax3.plot(ab_range, rrr_pred, 'g-', linewidth=2.5, label='Risk Reduction')
    ax3.axhline(y=50, color='orange', linestyle='--', linewidth=1.5, alpha=0.7, label='50% reduction')
    ax3.axhline(y=80, color='red', linestyle='--', linewidth=1.5, alpha=0.7, label='80% reduction')

    # Mark threshold for 50% and 80% reduction
    idx_50 = np.argmin(np.abs(rrr_pred - 50))
    idx_80 = np.argmin(np.abs(rrr_pred - 80))
    if rrr_pred[idx_50] >= 45:  # Only show if we achieve near 50%
        ax3.axvline(x=ab_range[idx_50], color='orange', linestyle=':', linewidth=1.5, alpha=0.7)
    if rrr_pred[idx_80] >= 75:  # Only show if we achieve near 80%
        ax3.axvline(x=ab_range[idx_80], color='red', linestyle=':', linewidth=1.5, alpha=0.7)

    ax3.set_xlabel(config.y_axis_label, fontsize=16, fontweight='bold')
    ax3.set_ylabel('Relative Risk Reduction (%)\n(vs minimum)', fontsize=16, fontweight='bold')
    ax3.tick_params(axis='both', labelsize=14)
    ax3.set_xlim(0, config.x_axis_max)
    ax3.set_ylim(0, 100)
    ax3.grid(True, alpha=0.3)
    ax3.legend(loc='lower right', frameon=True, fontsize=11)
    ax3.set_title('C. Relative Risk Reduction', fontsize=18, fontweight='bold', pad=10)

    plt.tight_layout()

    # Save with modified filename
    output_file = config.output_figure.replace('.pdf', '_OR_RRR.pdf')
    plt.savefig(output_file, dpi=FIGURE_DPI, bbox_inches='tight')
    plt.close()

    return output_file, reference_ab, or_pred, rrr_pred


def create_figure(config: AssayConfig, df_long: pd.DataFrame,
                  ab_range: np.ndarray, prob_pred: np.ndarray,
                  prob_lower: np.ndarray, prob_upper: np.ndarray):
    """Create publication-quality two-panel figure."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 10), gridspec_kw={'width_ratios': [0.6, 1]})
    
    # Panel A: Distribution
    np.random.seed(42)
    cases_data = df_long[df_long['outcome'] == 1]['antibody']
    controls_data = df_long[df_long['outcome'] == 0]['antibody']
    
    jitter_amount = 0.1
    x_controls = np.random.normal(0, jitter_amount, size=len(controls_data))
    x_cases = np.random.normal(1, jitter_amount, size=len(cases_data))
    
    ax1.scatter(x_controls, controls_data, alpha=0.6, s=40,
               color=CONTROL_COLOR, edgecolors='black', linewidth=0.5,
               label=f'Controls (n={len(controls_data)})')
    ax1.scatter(x_cases, cases_data, alpha=0.6, s=40,
               color=CASE_COLOR, edgecolors='black', linewidth=0.5,
               label=f'Cases (n={len(cases_data)})')
    
    ax1.hlines(stats.gmean(controls_data), -0.2, 0.2,
              colors='black', linewidth=3, label='Geometric mean')
    ax1.hlines(stats.gmean(cases_data), 0.8, 1.2,
              colors='black', linewidth=3)
    
    ax1.set_yscale('log')
    ax1.set_ylim(bottom=10**config.y_axis_bottom)
    ax1.set_xlim(-0.5, 1.5)
    ax1.set_xticks([0, 1])
    ax1.set_xticklabels(['Control', r'$\it{E. coli}$' + '\nsepsis'], fontsize=20, fontweight='bold')
    ax1.tick_params(axis='x', rotation=45)
    ax1.set_ylabel(config.y_axis_label, fontsize=20, fontweight='bold')
    ax1.tick_params(axis='y', labelsize=20)
    ax1.grid(True, alpha=0.3, axis='y')
    ax1.legend(loc='upper center', bbox_to_anchor=(0.5, -0.45), frameon=True, fontsize=15, ncol=2)
    
    # Panel B: Probability curve
    ax2.fill_between(ab_range, prob_lower, prob_upper,
                    alpha=0.5, color='gray', label='95% CI', zorder=1, edgecolor='none')
    ax2.plot(ab_range, prob_pred, 'k-', linewidth=2.5, label='Predicted probability', zorder=2)
    
    ax2.set_xlabel(config.y_axis_label, fontsize=24, fontweight='bold')
    ax2.set_ylabel('Probability of Sepsis', fontsize=24, fontweight='bold')
    ax2.tick_params(axis='both', labelsize=20)
    ax2.set_xlim(0, config.x_axis_max)
    ax2.set_ylim(0, max(prob_pred) * 1.1)
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc='upper right', frameon=True, fontsize=15)
    
    plt.tight_layout()
    plt.savefig(config.output_figure, dpi=FIGURE_DPI, bbox_inches='tight')
    plt.close()


def save_summary(config: AssayConfig, desc_stats: dict, matched_diffs: np.ndarray,
                 result, ab_range: np.ndarray, prob_pred: np.ndarray, 
                 baseline_risk: float):
    """Save analysis summary to text file."""
    coef = result.params['log_antibody']
    se = result.bse['log_antibody']
    ci_lower = result.conf_int().loc['log_antibody', 0]
    ci_upper = result.conf_int().loc['log_antibody', 1]
    p_value = result.pvalues['log_antibody']
    or_10fold = np.exp(coef)
    or_2fold = np.exp(coef * np.log10(2))
    
    with open(config.output_summary, 'w') as f:
        f.write(f"MATCHED CASE-CONTROL ANTIBODY ANALYSIS SUMMARY - {config.name}\n")
        f.write("=" * 70 + "\n\n")
        
        f.write("SAMPLE SIZE:\n")
        f.write(f"  Cases: {len(desc_stats['cases'])}\n")
        f.write(f"  Controls: {len(desc_stats['controls'])}\n")
        f.write(f"  Matched sets: {len(matched_diffs)}\n\n")
        
        f.write("DESCRIPTIVE STATISTICS:\n")
        f.write(f"  Cases - Geometric mean: {desc_stats['cases_gmean']:.1f}\n")
        f.write(f"  Controls - Geometric mean: {desc_stats['controls_gmean']:.1f}\n")
        f.write(f"  Geometric mean ratio: {desc_stats['gmr']:.3f}\n\n")
        
        f.write("CONDITIONAL LOGISTIC REGRESSION:\n")
        f.write(f"  Coefficient: {coef:.3f} (SE: {se:.3f})\n")
        f.write(f"  95% CI: [{ci_lower:.3f}, {ci_upper:.3f}]\n")
        f.write(f"  P-value: {p_value:.4e}\n\n")
        
        f.write("ODDS RATIOS:\n")
        f.write(f"  Per 10-fold increase: {or_10fold:.3f}\n")
        f.write(f"  Per 2-fold increase: {or_2fold:.3f}\n\n")
        
        f.write("PROTECTIVE THRESHOLDS:\n")
        for reduction in [0.50, 0.75, 0.90, 0.95]:
            target_risk = baseline_risk * (1 - reduction)
            idx = np.argmin(np.abs(prob_pred - target_risk))
            threshold = ab_range[idx]
            f.write(f"  {reduction*100:.0f}% risk reduction: {threshold:.1f}\n")


def run_analysis(config: AssayConfig, verbose: bool = True):
    """Run complete analysis for a single assay."""
    if verbose:
        print("\n" + "=" * 70)
        print(f"MATCHED CASE-CONTROL ANTIBODY ANALYSIS: {config.name}")
        print("=" * 70)
    
    # Load data
    df_long = load_and_preprocess_data(config.data_str)
    if verbose:
        print(f"\nLoaded {df_long['match_id'].nunique()} matched sets")
        print(f"  Total observations: {len(df_long)}")
        print(f"  Cases: {df_long['outcome'].sum()}")
        print(f"  Controls: {len(df_long) - df_long['outcome'].sum()}")
    
    # Descriptive statistics
    desc_stats = compute_descriptive_stats(df_long)
    matched_diffs = compute_matched_differences(df_long)
    
    if verbose:
        print(f"\nDescriptive Statistics:")
        print(f"  Cases - Geometric mean: {desc_stats['cases_gmean']:.1f}")
        print(f"  Controls - Geometric mean: {desc_stats['controls_gmean']:.1f}")
        print(f"  Geometric mean ratio: {desc_stats['gmr']:.3f}")
    
    # Conditional logistic regression
    result = fit_conditional_logistic(df_long)
    
    coef = result.params['log_antibody']
    se = result.bse['log_antibody']
    ci_lower = result.conf_int().loc['log_antibody', 0]
    ci_upper = result.conf_int().loc['log_antibody', 1]
    p_value = result.pvalues['log_antibody']
    or_10fold = np.exp(coef)
    
    if verbose:
        print(f"\nConditional Logistic Regression:")
        print(f"  Coefficient: {coef:.3f} (SE: {se:.3f})")
        print(f"  95% CI: [{ci_lower:.3f}, {ci_upper:.3f}]")
        print(f"  P-value: {p_value:.4e}")
        print(f"  OR per 10-fold increase: {or_10fold:.3f}")
    
    # Generate predictions
    if verbose:
        print("\nCalculating bootstrap confidence intervals...")
    
    if config.use_standard_logistic_for_viz:
        ab_range, prob_pred, prob_lower, prob_upper, successful_boots = \
            generate_predictions_standard_logistic(df_long)
    else:
        ab_range, prob_pred, prob_lower, prob_upper, successful_boots = \
            generate_predictions_conditional(df_long, result)
    
    if verbose:
        print(f"  Bootstrap success rate: {successful_boots}/1000")
    
    # Create figures
    if verbose:
        print("\nGenerating figures...")
    create_figure(config, df_long, ab_range, prob_pred, prob_lower, prob_upper)

    # Create OR/RRR figure (more interpretable for case-control)
    if verbose:
        print("Generating Odds Ratio / Risk Reduction figure...")
    or_rrr_file, reference_ab, or_pred, rrr_pred = create_figure_or_rrr(
        config, df_long, ab_range, result
    )

    # Save summary
    baseline_risk = df_long['outcome'].mean()
    save_summary(config, desc_stats, matched_diffs, result, ab_range, prob_pred, baseline_risk)

    if verbose:
        print(f"\nOutputs generated:")
        print(f"  - {config.output_figure}")
        print(f"  - {or_rrr_file}")
        print(f"  - {config.output_summary}")
        print(f"\nOR/RRR interpretation (reference = minimum: {reference_ab:.1f}):")
        # Find key thresholds
        idx_50 = np.argmin(np.abs(rrr_pred - 50))
        idx_80 = np.argmin(np.abs(rrr_pred - 80))
        print(f"  50% risk reduction at: {ab_range[idx_50]:.1f}")
        print(f"  80% risk reduction at: {ab_range[idx_80]:.1f}")
    
    return {
        'config': config,
        'df_long': df_long,
        'result': result,
        'desc_stats': desc_stats,
        'ab_range': ab_range,
        'prob_pred': prob_pred,
    }


# =============================================================================
# MAIN ENTRY POINT
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Matched Case-Control Analysis of Maternal Antibodies',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python %(prog)s                    # Run all assays
    python %(prog)s --assay THP1       # Run single assay
    python %(prog)s --assay Mix8 HL60  # Run multiple assays
    python %(prog)s --list             # List available assays
        """
    )
    parser.add_argument(
        '--assay', '-a',
        nargs='+',
        choices=list(ASSAY_CONFIGS.keys()),
        help='Assay(s) to analyze (default: all)'
    )
    parser.add_argument(
        '--list', '-l',
        action='store_true',
        help='List available assays and exit'
    )
    parser.add_argument(
        '--quiet', '-q',
        action='store_true',
        help='Suppress verbose output'
    )
    
    args = parser.parse_args()
    
    if args.list:
        print("\nAvailable assays:")
        for key, config in ASSAY_CONFIGS.items():
            print(f"  {key:6} - {config.name} ({config.y_axis_label})")
        return
    
    setup_plot_style()
    
    assays_to_run = args.assay if args.assay else list(ASSAY_CONFIGS.keys())
    
    for assay_key in assays_to_run:
        config = ASSAY_CONFIGS[assay_key]
        run_analysis(config, verbose=not args.quiet)
    
    print("\n" + "=" * 70)
    print("ALL ANALYSES COMPLETE")
    print("=" * 70)


if __name__ == '__main__':
    main()
