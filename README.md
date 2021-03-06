# parallel-programming-3

## Requirements
* CMake >= 3.12.4

## Server/Master
サーバ/マスターの役割を担うプログラム（以下マスタープログラムと呼ぶ）。
クライアントへ問題の要求や解の送信を行う。
またワーカへの問題配布やワーカからの解の受信する。

クライアントプログラムを起動し、
マスタープログラム実行後、クライアントへの接続を行い、
指定されたワーカ数分、接続があるまで待機し、接続完了後問題の要求を開始する。
(クライアントプログラムは直接連絡をいただければ共有します)

### How to compile
```
cd coop/master 
cmake -DCMAKE_BUILD_TYPE=Release
make
```
### How to execution
マスタープログラムは以下のように指定する。
実行すると各行列乗算の実行時間が出力される。
```
./master クライアントのIPアドレス クライアントのポート番号 ワーカの向けポート番号 ワーカ数
```

## Worker
ワーカプログラム。
マスタープログラムから問題の受信や行列乗算、LU分解、解計算を行う。
解計算後、マスタープログラムへ解を送信する。

### How to compile
```
cd coop/worker
cmake -DCMAKE_BUILD_TYPE=Release
make
```
### How to execution
```
./worker サーバのIPアドレス サーバのポート番号
```

## MatMul
SIMD命令を使った行列乗算の実装
* 単純な行列乗算
* ブロック化した行列乗算
* ブロック化しスレッドを使った行列乗算
* シュトラッセンのアルゴリズムを使った行列乗算

シュトラッセンのアルゴリズムを使った行列乗算では、単純な行列乗算と比較すると誤差が出ている。
シュトラッセンのアルゴリズムでは行列の減算を行っており、桁落ちが発生していると考えている。

[シュトラッセンのアルゴリズム](https://ja.wikipedia.org/wiki/%E3%82%B7%E3%83%A5%E3%83%88%E3%83%A9%E3%83%83%E3%82%BB%E3%83%B3%E3%81%AE%E3%82%A2%E3%83%AB%E3%82%B4%E3%83%AA%E3%82%BA%E3%83%A0)

### How to compile
```
cd matmul
cmake -DCMAKE_BUILD_TYPE=Release
make
```

### How to execution
行列のサイズは以下のように指定する。実行すると各行列乗算の実行時間が出力される。
```
./matmul 行列のサイズ
```

## LU分解
ブロック形式ガウス方のLU分解の実装。
実装が間違っていたので一時的に削除。

P.33
[参考PDF](http://www.cspp.cc.u-tokyo.ac.jp/hanawa/class/spc2016s/sp20160614-2.pdf)
