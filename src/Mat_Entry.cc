#include "Mat_Entry.h"

Entry operator+ (const Entry& ent_1, const Entry& ent_2) {
    assert(ent_1.getSize()==ent_2.getSize());
    Entry ent_ot(ent_1.getSize());
    ent_ot.getVal() = ent_1.getVal() + ent_2.getVal();
    return ent_ot;
}

Entry operator+ (const Entry& ent_1, const DType& val) {
    Entry ent_ot(ent_1.getSize());
    ent_ot.getVal() = ent_1.getVal() + val;
    return ent_ot;
}

Entry operator+ (const DType& val, const Entry& ent_2) {
    Entry ent_ot(ent_2.getSize());
    ent_ot.getVal() = val + ent_2.getVal();
    return ent_ot;
}

Entry operator+ (const Entry& ent_1) {
    Entry ent_ot(ent_1.getSize());
    ent_ot.getVal() = ent_1.getVal();
    return ent_ot;
}

Entry operator- (const Entry& ent_1, const Entry& ent_2) {
    assert(ent_1.getSize()==ent_2.getSize());
    Entry ent_ot(ent_1.getSize());
    ent_ot.getVal() = ent_1.getVal() - ent_2.getVal();
    return ent_ot;
}

Entry operator- (const Entry& ent_1, DType val) {
    Entry ent_ot(ent_1.getSize());
    ent_ot.getVal() = ent_1.getVal() - val;
    return ent_ot;
}

Entry operator- (DType val, const Entry& ent_2) {
    Entry ent_ot(ent_2.getSize());
    ent_ot.getVal() = val - ent_2.getVal();
    return ent_ot;
}

Entry operator- (const Entry& ent_1) {
    Entry ent_ot(ent_1.getSize());
    ent_ot.getVal() = -ent_1.getVal();
    return ent_ot;
}

Entry operator* (const Entry& ent_1, const Entry& ent_2) {
    assert(ent_1.getSize()==ent_2.getSize());
    Entry ent_ot(ent_1.getSize());
    ent_ot.getVal() = ent_1.getVal() * ent_2.getVal();
    return ent_ot;
}

Entry operator* (const DType& sca_1, const Entry& ent_2) {
    Entry ent_ot(ent_2.getSize());
    ent_ot.getVal() = sca_1 * ent_2.getVal();
    return ent_ot;
}

Entry operator* (const Entry& ent_1, const DType& sca_2) {
    Entry ent_ot(ent_1.getSize());
    ent_ot.getVal() = sca_2 * ent_1.getVal();
    return ent_ot;
}

Entry operator/ (const Entry& ent_1, const DType& sca_2) {
    Entry ent_ot(ent_1.getSize());
    ent_ot.getVal() = ent_1.getVal()/sca_2;
    return ent_ot;
}

VComp operator+ (const VComp& comp_1, const VComp& comp_2) {
    assert(comp_1.getSize()==comp_2.getSize());
    VComp comp_ot(comp_1.getSize());
    comp_ot.getVal() = comp_1.getVal() + comp_2.getVal();
    return comp_ot;
}

VComp operator+ (const VComp& comp_1, const DType& val) {
    VComp comp_ot(comp_1.getSize());
    comp_ot.getVal() = comp_1.getVal() + val;
    return comp_ot;
}

VComp operator+ (const DType& val, const VComp& comp_2) {
    VComp comp_ot(comp_2.getSize());
    comp_ot.getVal() = val + comp_2.getVal();
    return comp_ot;
}

VComp operator+ (const VComp& comp_1) {
    VComp comp_ot(comp_1.getSize());
    comp_ot.getVal() = comp_1.getVal();
    return comp_ot;
}

VComp operator- (const VComp& comp_1, const VComp& comp_2) {
    assert(comp_1.getSize()==comp_2.getSize());
    VComp comp_ot(comp_1.getSize());
    comp_ot.getVal() = comp_1.getVal() - comp_2.getVal();
    return comp_ot;
}

VComp operator- (const VComp& comp_1, const DType& val) {
    VComp comp_ot(comp_1.getSize());
    comp_ot.getVal() = comp_1.getVal() - val;
    return comp_ot;
}

VComp operator- (const DType& val, const VComp& comp_2) {
    VComp comp_ot(comp_2.getSize());
    comp_ot.getVal() = val - comp_2.getVal();
    return comp_ot;
}

VComp operator- (const VComp& comp_1) {
    VComp comp_ot(comp_1.getSize());
    comp_ot.getVal() = -comp_1.getVal();
    return comp_ot;
}

// Inner product
DType operator* (const VComp& comp_1, const VComp& comp_2) {
    return comp_1.getVal() * comp_2.getVal();
}

VComp operator* (const VComp& comp_1, const DType& val) {
    VComp comp_ot(comp_1.getSize());
    comp_ot.getVal() = comp_1.getVal() * val;
    return comp_ot;
}

VComp operator* (const DType& val, const VComp& comp_2) {
    VComp comp_ot(comp_2.getSize());
    comp_ot.getVal() = val * comp_2.getVal();
    return comp_ot;
}

VComp operator* (const Entry& ent, const VComp& comp) {
    assert(ent.getSize()==comp.getSize());
    VComp comp_ot(comp.getSize());
    comp_ot.getVal() = ent.getVal() * comp.getVal();
    return comp_ot;
}

VComp operator/ (const VComp& comp_1, const DType& val) {
    VComp comp_ot(comp_1.getSize());
    comp_ot.getVal() = comp_1.getVal() / val;
    return comp_ot;
}
