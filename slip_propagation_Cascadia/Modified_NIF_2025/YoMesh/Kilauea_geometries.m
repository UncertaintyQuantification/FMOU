[origin] = [-155.273256184167 19.339057464198 ]; % Site of MANE

east_rift = [ 700 6815;
    2160 5697
    3001 4989
    4198 4148 
    5512 3429
    6646 2842
    8664 2799
    10230 2744
    11670 2834
    13560 3134
    17500 5300
    20320 6700
    21980 7093
    24340 8436
    26360 9568
    28900 10340
    29860 10800
    33110 11160
    37020 12830
    40300 15480
    45000 18000];
east_rift = curvspace(east_rift,round((east_rift(end,1)-east_rift(1,1))/1000));    

SW_rift = [-2231 7081
    -3418 6107
    -5908 3650
    -8557 784
    -11280 -1295
    -13230 -5488
    -14330 -7734
    -15350 -10720
    -15910 -13930
    -16640 -16550
    -17150 -18510];
SW_rift = curvspace(SW_rift,abs(round((SW_rift(end,1)-SW_rift(1,1))/1000)));    

hilina = [-8532 -10870
    -7454 -8171
    -6326 -7308
    -4734 -5763
    -3206 -4940
    -1384 -4390
    61 -3781
    2549 -3054
    5829 -2748
    6850 -2394
    9563 -2552
    12070 -2692];
hilina = curvspace(hilina,round((hilina(end,1)-hilina(1,1))/1000));

holei = [3000 -7750
    3948 -7000
    6323 -6190
    8038 -4918 
    9621 -3939
    10500 -3400
    11540 -3040
    13000 -3000
    14070 -2987
    16000 -3012
    17850 -2493
    20420 -2112
    22680 -1738];
holei = curvspace(holei,round((holei(end,1)-holei(1,1))/1000));

koae = [7036 1758
    3117 1248
    1141 893
    -308 306
    -3661 -600
    -6015 -1408 
    -7171 -2083
    -9589 -3476];
koae = curvspace(koae,abs(round((koae(end,1)-koae(1,1))/1000)));

plot(line_xy(1,:),line_xy(2,:),'k'); hold on; 
axis equal; axis([-20 50 -20 20]*1000); 
plot(SW_rift(:,1),SW_rift(:,2),'go-','Linewidth',2)
plot(east_rift(:,1),east_rift(:,2),'bo-','Linewidth',2)
plot(hilina(:,1), hilina(:,2),'ro-','Linewidth',2)
plot(holei(:,1), holei(:,2),'co-','Linewidth',2)
plot(koae(:,1), koae(:,2),'mo-','Linewidth',2)
title('origin = [-155.273256184167 19.339057464198]')
xlabel('East (km)'); ylabel('North (km)')