#include <stdio.h>
#include <iostream>
#include <memory>
#include <cmath>

// �v�Z�O���t�̃m�[�h�N���X---------------------------------------------------------------
class Node {
public:

    // �f�[�^�����o-----------------------------------------------------------------

    double v;   // �֐��l
    double d;   // �����l

    // �q�m�[�h�ւ̃|�C���^
    std::shared_ptr<Node> L_Node;
    std::shared_ptr<Node> R_Node;

    // �q�m�[�h�̗v�f�I�Γ��֐��l
    double L_d;
    double R_d;

    // �R���X�g���N�^---------------------------------------------------------------
    Node(double val,                               // �֐��l
        std::shared_ptr<Node> l_node = nullptr,    // �q�m�[�h�ւ̃|�C���^
        std::shared_ptr<Node> r_node = nullptr,    // �q�m�[�h�ւ̃|�C���^
        double l_d = 0.0,                          // �q�m�[�h�ւ̗v�f�I�Γ��֐��l
        double r_d = 0.0)                          // �q�m�[�h�ւ̗v�f�I�Γ��֐��l
        : v(val), d(0.0), L_Node(l_node), R_Node(r_node), L_d(l_d), R_d(r_d) {}

    // �����l�̌v�Z-----------------------------------------------------------------
    void backward(double parent_d, double parent_elemental_d) {

        // ���g���ŉ��w���`�F�b�N
        if (L_Node != nullptr) {  

            // �ŉ��w�łȂ���Δ����l���X�V����
            d = parent_d * parent_elemental_d;

            // �q�m�[�h�ɑ΂��ċA�Ăяo�����s��
            L_Node->backward(d, L_d);
            if (R_Node != nullptr) R_Node->backward(d, R_d);

        // ���g���ŉ��w�̏ꍇ
        } else {                  

            // �����l�𑫂�����
            d += parent_d * parent_elemental_d;
        }
    }
};

// Node����p�̃N���X---------------------------------------------------------------------
class Var {
public:

    // �f�[�^�����o-----------------------------------------------------------------
    std::shared_ptr<Node> p; // �m�[�h�ւ̃|�C���^

    // �R���X�g���N�^---------------------------------------------------------------
    Var(double value = 0.0, 
        std::shared_ptr<Node> node_L = nullptr, 
        std::shared_ptr<Node> node_R = nullptr, 
        double ele_d_L = 0.0, 
        double ele_d_R = 0.0)
        : p(std::make_shared<Node>(value, node_L, node_R, ele_d_L, ele_d_R)) {}

    // �֐��l�̎擾
    double value() const {
        return p->v;
    }
    // �����l�̎擾
    double derivative() const {
        return p->d;
    }
    // �����l�̌v�Z
    void backward() {
        p->d = 1.0;
        p->backward(1.0, 1.0);
    }

    // ���Z�q�I�[�o�[���[�h---------------------------------------------------------

    // ���Z�@�i�������e�̐����̂��߂ɏڍׂɓ��e���L�q�j
    friend Var operator+(const Var& a, const Var& b) {

        Var r;                            // ���Z���ʂ̃m�[�h�𐶐�
        r.p->v = a.value() + b.value();   // ���Z�����l���Z�b�g
        r.p->d = 0.0;                     // ���Z���ʃm�[�h�̔����l��������
        r.p->L_Node = a.p;                // �q�m�[�h(��)�ւ̃|�C���^���Z�b�g
        r.p->R_Node = b.p;                // �q�m�[�h(�E)�ւ̃|�C���^���Z�b�g
        r.p->L_d = 1.0;                   // ���Z�ɂ��v�f�I�Γ��֐��l���Z�b�g
        r.p->R_d = 1.0;                   // ���Z�ɂ��v�f�I�Γ��֐��l���Z�b�g

        return r;
    }
    // ���Z�@�i�ȉ��ł́C��L���Z�����Ɠ��l�̓��e��1�s�ŋL�q�j
    friend Var operator-(const Var& a, const Var& b) {
        return Var(a.value() - b.value(), a.p, b.p, 1.0, -1.0);
    }
    // ��Z
    friend Var operator*(const Var& a, const Var& b) {
        return Var(a.value() * b.value(), a.p, b.p, b.value(), a.value());
    }
    // ���Z 
    friend Var operator/(const Var& a, const Var& b) {
        return Var(a.value() / b.value(), a.p, b.p, 
            1.0 / b.value(), - a.value() / (b.value()*b.value()));
    }
    // �P���}�C�i�X���Z�q
    friend Var operator-(const Var& a) {
        return Var(a.value(), a.p, nullptr, -1.0, 0.0);
    }

    // �����֐�
    friend Var sin(const Var& a) {
        return Var(std::sin(a.value()), a.p, nullptr, std::cos(a.value()), 0.0);
    }
    friend Var exp(const Var& a) {
        return Var(std::exp(a.value()), a.p, nullptr, std::exp(a.value()), 0.0);
    }
};

// ���C���֐�-----------------------------------------------------------------------------
int main(){

    Var x(2.0);

    Var f = sin(x * x * x - exp(x)) / (x * x + 2.0);

    f.backward();

    printf("f(%.1lf) = %lf\n", x.value(), f.value());
    printf("f'(%.1lf) = %lf\n", x.value(), x.derivative());
    
    return 0;
}
