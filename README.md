# a-analog_MZV
culculate a-analog muliple zeta value
https://yuttarikanata.github.io/a-analog_MZV/

まず、やり方を残しておくと後世のためになるよ
emccでgmpをアッセンブリにするときにちゃんと
「コンパイラがgccではなくemccになっているか」
を確認しよう
(emconfigure ./configure --host=wasm32-unknown-emscripten --prefix=${HOME}/optのやつ)
参考にしたpdf: https://www.kurims.kyoto-u.ac.jp/~kyodo/kokyuroku/contents/pdf/2273-08.pdf

そして、それができない場合はディレクトリを新しく作ったやつにしよう
実際にやったのは/optを放棄して/gmp-bulidをつくってその下に/gmp-build/gmp-6.3.0みたいに入れてやってた
これしたらpermission deniedみたいなやつが消えた

むやみにフォルダを消すのはやめよう
rmコマンド強制的に使ってoptファイルを消したらなぜかいろんなコマンドが使えなくなっていた！(curlとか)
WSL2入れなおす羽目になった

ちゃんとコードは見よう
特にメモリのところ。
GMPのmpf_t型の配列作ったところとか変数作ったところとかちゃんと開放してるかとか
あとはemccでwasmにするときに再帰を深くする場合はSTACK_SIZE=128MBとかちゃんとスタックサイズふやそう
これだけのせいで2日ぐらい消えた。index out of boundsって言われるからコードに問題があると思うやん！なかった

以下Chat GPT


あなたのコンパイルコマンドは、ほとんどの場合、問題なく動作するはずですが、いくつかのポイントを確認しておくとより万全です。
コマンドの詳細な確認

    bash
    
    emcc Azeta.c -o Azeta.js \
      -I${HOME}/gmp-build/include \
      -static \
      -L${HOME}/gmp-build/lib \
      -lgmp \
      -s WASM=1 \
      -s EXPORTED_FUNCTIONS='["_azeta","_malloc","_free"]' \
      -s EXPORTED_RUNTIME_METHODS='["ccall","cwrap","setValue","getValue"]' \
      -s INITIAL_MEMORY=1024MB \
      -s ALLOW_MEMORY_GROWTH=1 \
      -lm \
      -O3 \
      -s STACK_SIZE=1MB

チェックポイント

    -static の使用:
        -staticオプションは、静的リンクを強制しますが、Emscriptenの標準的な使用ではあまり必要ありません。
        EmscriptenはWebAssemblyのために標準で動的リンクライクに扱うので、-staticは省略しても大丈夫です。

    **INITIAL_MEMORY=1024MBとALLOW_MEMORY_GROWTH=1**の組み合わせ:
        INITIAL_MEMORY=1024MBは大きなメモリを最初に確保し、
        その後ALLOW_MEMORY_GROWTH=1によりさらにメモリを増やすことが可能です。
        この設定は良いですが、初期メモリサイズが大きいと、プログラムが最初から多くのメモリを要求し、
        ブラウザでのロードが遅くなることがあるため注意が必要です。必要ならばメモリサイズを少し下げて、
        成長オプションに依存するようにしても良いでしょう。

    STACK_SIZE=1MB:
        1MBのスタックサイズは通常のアプリケーションには十分ですが、大きな局所変数や深い再帰呼び出しがある場合は、
        スタックオーバーフローのリスクがあります。そのような場合、スタックサイズを増やすことを検討してください。

    -lm の位置:
        -lmオプションはリンクオプションの後に置く必要がありますが、既に正しい位置にあるので問題ありません。

    最適化オプション -O3:
        -O3は最大の最適化を行いますが、場合によってはコードサイズが大きくなることがあります。
        もしサイズが問題になる場合は、-Os（サイズ最適化）を試すことも検討できます。

その他の考慮点

    ブラウザ互換性: INITIAL_MEMORYが大きい場合、ブラウザのメモリ制限に引っかかる可能性があるため、複数のブラウザでテストして互換性を確認すると良いでしょう。

    デバッグ: もし動作に問題がある場合、-gオプションを追加してデバッグ情報を含め、EmscriptenのASSERTIONS=1オプションを使って問題を特定しやすくすることができます。

全体的には、このコマンドは高いパフォーマンスを意識した設定となっており、万全に近いと言えますが、上記の点を考慮に入れてさらに調整できる部分があるかもしれません。
