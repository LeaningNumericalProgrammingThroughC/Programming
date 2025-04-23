#include <stdio.h>
#include <iostream>
#include <memory>
#include <cmath>

// 計算グラフのノードクラス---------------------------------------------------------------
class Node {
public:

    // データメンバ-----------------------------------------------------------------

    double v;   // 関数値
    double d;   // 微分値

    // 子ノードへのポインタ
    std::shared_ptr<Node> L_Node;
    std::shared_ptr<Node> R_Node;

    // 子ノードの要素的偏導関数値
    double L_d;
    double R_d;

    // コンストラクタ---------------------------------------------------------------
    Node(double val,                               // 関数値
        std::shared_ptr<Node> l_node = nullptr,    // 子ノードへのポインタ
        std::shared_ptr<Node> r_node = nullptr,    // 子ノードへのポインタ
        double l_d = 0.0,                          // 子ノードへの要素的偏導関数値
        double r_d = 0.0)                          // 子ノードへの要素的偏導関数値
        : v(val), d(0.0), L_Node(l_node), R_Node(r_node), L_d(l_d), R_d(r_d) {}

    // 微分値の計算-----------------------------------------------------------------
    void backward(double parent_d, double parent_elemental_d) {

        // 自身が最下層かチェック
        if (L_Node != nullptr) {  

            // 最下層でなければ微分値を更新する
            d = parent_d * parent_elemental_d;

            // 子ノードに対し再帰呼び出しを行う
            L_Node->backward(d, L_d);
            if (R_Node != nullptr) R_Node->backward(d, R_d);

        // 自身が最下層の場合
        } else {                  

            // 微分値を足し込む
            d += parent_d * parent_elemental_d;
        }
    }
};

// Node操作用のクラス---------------------------------------------------------------------
class Var {
public:

    // データメンバ-----------------------------------------------------------------
    std::shared_ptr<Node> p; // ノードへのポインタ

    // コンストラクタ---------------------------------------------------------------
    Var(double value = 0.0, 
        std::shared_ptr<Node> node_L = nullptr, 
        std::shared_ptr<Node> node_R = nullptr, 
        double ele_d_L = 0.0, 
        double ele_d_R = 0.0)
        : p(std::make_shared<Node>(value, node_L, node_R, ele_d_L, ele_d_R)) {}

    // 関数値の取得
    double value() const {
        return p->v;
    }
    // 微分値の取得
    double derivative() const {
        return p->d;
    }
    // 微分値の計算
    void backward() {
        p->d = 1.0;
        p->backward(1.0, 1.0);
    }

    // 演算子オーバーロード---------------------------------------------------------

    // 加算　（処理内容の説明のために詳細に内容を記述）
    friend Var operator+(const Var& a, const Var& b) {

        Var r;                            // 加算結果のノードを生成
        r.p->v = a.value() + b.value();   // 加算した値をセット
        r.p->d = 0.0;                     // 加算結果ノードの微分値を初期化
        r.p->L_Node = a.p;                // 子ノード(左)へのポインタをセット
        r.p->R_Node = b.p;                // 子ノード(右)へのポインタをセット
        r.p->L_d = 1.0;                   // 加算による要素的偏導関数値をセット
        r.p->R_d = 1.0;                   // 加算による要素的偏導関数値をセット

        return r;
    }
    // 減算　（以下では，上記加算処理と同様の内容を1行で記述）
    friend Var operator-(const Var& a, const Var& b) {
        return Var(a.value() - b.value(), a.p, b.p, 1.0, -1.0);
    }
    // 乗算
    friend Var operator*(const Var& a, const Var& b) {
        return Var(a.value() * b.value(), a.p, b.p, b.value(), a.value());
    }
    // 除算 
    friend Var operator/(const Var& a, const Var& b) {
        return Var(a.value() / b.value(), a.p, b.p, 
            1.0 / b.value(), - a.value() / (b.value()*b.value()));
    }
    // 単項マイナス演算子
    friend Var operator-(const Var& a) {
        return Var(a.value(), a.p, nullptr, -1.0, 0.0);
    }

    // 初等関数
    friend Var sin(const Var& a) {
        return Var(std::sin(a.value()), a.p, nullptr, std::cos(a.value()), 0.0);
    }
    friend Var exp(const Var& a) {
        return Var(std::exp(a.value()), a.p, nullptr, std::exp(a.value()), 0.0);
    }
};

// メイン関数-----------------------------------------------------------------------------
int main(){

    Var x(2.0);

    Var f = sin(x * x * x - exp(x)) / (x * x + 2.0);

    f.backward();

    printf("f(%.1lf) = %lf\n", x.value(), f.value());
    printf("f'(%.1lf) = %lf\n", x.value(), x.derivative());
    
    return 0;
}
