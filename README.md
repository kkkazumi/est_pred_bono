# mental estimation prog tutorial

#### 熊谷修論実験
#### クイズゲームで得られたデータから気分推定を行うプログラム

## 内容：
#### README.md
* 説明が書いてあります。

#### factor.csv
* 因子群のデータ

	+ クイズゲームの状況とロボット行動を表す

	+ 各変数
		+ クイズ回数、正答率、共感行動実施率、励まし行動実施率、からかい行動実施率、クイズと関係ない行動実施率、発話・動作無し確率、獲得点数、連勝数、連敗数

#### signal.csv

* 表情認識結果のデータ群

	- 喜び、驚き、怒り、悲しみの表情が全体の表情に占める割合を表す
	- 各変数
		- 喜び、驚き、怒り、悲しみ

#### estimate.cpp

* `factor.csv`と`signal.csv`から気分推定結果の時系列データを算出するプログラム
* 算出した気分推定結果の時系列データは`mental.csv` に出力される

#### predict.cpp
* `factor.csv`と`mental.csv`から気分変化予測結果の時系列データを算出するプログラム
* 算出した気分変化予測結果の時系列データはあるファイル(名前はまだ無い)に出力される

## プログラムを動かす準備
* Eigenをインストール

	```	
	sudo apt-get install libeigen3-dev #バージョンは任意
	```

*  gccの最新版をインストールしておく

	```
	$sudo apt-get install g++-4.8
	```

## プログラムのうごかし方

1. `estimate.cpp`のコンパイル

	```
	$g++ estimate.cpp -std=c++11 -o estimate
	```

2. `estimate` を実行する
	- 気分推定結果を出力するファイル名を引数にする

	```
	$./estimate test
	```
	
	- `/test

## 気分推定結果と被験者自身による気分評価結果との相関計算

### sokan.R

1. `R` を立ち上げる

2. `sokan.R` を読み込む。 `source('sokan.R')`

3. `q()` で `R` からシェルに戻る。

4. `soukan.csv` 内に出力されているよ。#(解説必要なくね・・？)


## func_new.py

- contains new functions to calculate again.#mental estimation prog tutorial

## face vs mood
- as the master thesis project, two methods of selecting robot's behavior based on 1)face and 2)mood are compared.
  - note: orders of counter balance
    - 1: 1)mood,2)face
    - 2: 1)mood,2)face
    - 3: 1)mood, 2)face
    - 4: 1)face,2)mood
    - 5: 1)face, 2)mood
    - 6: 1)mood, 2)face
    - 7: 1)mood, 2)face
    - 8: 1)mood, 2)face
    - 9: 1)face,2)mood
