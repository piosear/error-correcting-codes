# Turbo codes
通信規格LTEなどに用いられている誤り訂正符号であるTurbo符号の復号誤り率を求めるシミュレーションを行うプログラムです。
## Description
このシミュレーションではベーシックなTurbo符号を用いています。  
変調方式としてBPSK (binary phase-shift keying)、復号法としてLog-MAP復号を用いています。  
生成されるファイルはビット誤り率(BER: bit error rate)特性とフレーム誤り率(FER: frame error rate)特性の2つです。 
変更可能なパラメータは  
- INFO_LENGTH
- NU (2,3,4から選択可)  
- ITERATE  

それぞれ増やすほどに誤り率特性が向上しますが、計算量が増大します．  

情報長1000、レジスタ数3、反復回数8を基本として、１つのパラメータを変化させたときのBER特性のグラフは次のようになります。
<img src="https://github.com/piosear/error-correcting-codes/blob/master/turbo-codes/Figures/Turbo_performance.png" alt="graph" title="Turbo符号の誤り率特性" width=50%>

横軸はビット電力対雑音電力比<img src="https://latex.codecogs.com/gif.latex?E_b/N_0" />
を表し、右にいくほど大きな電力で送信していることを示します。  
つまり、左下にいくほど小さな電力で小さな誤り率を達成しているので、特性が良いことになります。  


