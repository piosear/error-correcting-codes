# Turbo codes
このプログラムは、通信規格LTEなどに用いられている誤り訂正符号であるTurbo符号の復号誤り率を求めるシミュレーションです。
## Description
このプログラムはベーシックなTurbo符号のシミュレーションです。  
変調方式としてBPSK (binary phase-shift keying)、復号法としてLog-MAP復号を用いています。  
生成されるファイルは横軸がビット電力対雑音電力比$E_b/N_0$、縦軸がビット誤り率(BER: bit error rate)のBER特性と、横軸が$E_b/N_0$で縦軸がフレーム誤り率(FER: frame error rate)のFER特性の2つです。  
変更可能なパラメータは  
- INFO_LENGTH
- NU (2,3,4から選択可)  
- ITERATE  

それぞれ増やすほどに誤り率特性が向上しますが、計算量が増大します．  

情報長1000、レジスタ数3、反復回数8を基本として、１つのパラメータを変化させたときのグラフは次のようになります。

<img src="https://github.com/piosear/error-correcting-codes/blob/master/turbo-codes/Figures/Turbo_performance.png" alt="graph" title="Turbo符号の誤り率特性">
