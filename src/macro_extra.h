#define DESIG_1(v) [0]=v
#define DESIG__1(n,v) [0]={DESIG_##n(v)}
#define EmptyPyGetSetDef_1 [0]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_1 [0]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_1 [0]={NULL,0,0,0,NULL}
#define DESIG_2(v) DESIG_1(v),[1]=v
#define DESIG__2(n,v) DESIG__1(n,v),[1]={DESIG_##n(v)}
#define EmptyPyGetSetDef_2 EmptyPyGetSetDef_1,[1]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_2 EmptyPyMethodDef_1,[1]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_2 EmptyPyMemberDef_1,[1]={NULL,0,0,0,NULL}
#define DESIG_3(v) DESIG_2(v),[2]=v
#define DESIG__3(n,v) DESIG__2(n,v),[2]={DESIG_##n(v)}
#define EmptyPyGetSetDef_3 EmptyPyGetSetDef_2,[2]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_3 EmptyPyMethodDef_2,[2]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_3 EmptyPyMemberDef_2,[2]={NULL,0,0,0,NULL}
#define DESIG_4(v) DESIG_3(v),[3]=v
#define DESIG__4(n,v) DESIG__3(n,v),[3]={DESIG_##n(v)}
#define EmptyPyGetSetDef_4 EmptyPyGetSetDef_3,[3]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_4 EmptyPyMethodDef_3,[3]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_4 EmptyPyMemberDef_3,[3]={NULL,0,0,0,NULL}
#define DESIG_5(v) DESIG_4(v),[4]=v
#define DESIG__5(n,v) DESIG__4(n,v),[4]={DESIG_##n(v)}
#define EmptyPyGetSetDef_5 EmptyPyGetSetDef_4,[4]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_5 EmptyPyMethodDef_4,[4]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_5 EmptyPyMemberDef_4,[4]={NULL,0,0,0,NULL}
#define DESIG_6(v) DESIG_5(v),[5]=v
#define DESIG__6(n,v) DESIG__5(n,v),[5]={DESIG_##n(v)}
#define EmptyPyGetSetDef_6 EmptyPyGetSetDef_5,[5]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_6 EmptyPyMethodDef_5,[5]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_6 EmptyPyMemberDef_5,[5]={NULL,0,0,0,NULL}
#define DESIG_7(v) DESIG_6(v),[6]=v
#define DESIG__7(n,v) DESIG__6(n,v),[6]={DESIG_##n(v)}
#define EmptyPyGetSetDef_7 EmptyPyGetSetDef_6,[6]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_7 EmptyPyMethodDef_6,[6]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_7 EmptyPyMemberDef_6,[6]={NULL,0,0,0,NULL}
#define DESIG_8(v) DESIG_7(v),[7]=v
#define DESIG__8(n,v) DESIG__7(n,v),[7]={DESIG_##n(v)}
#define EmptyPyGetSetDef_8 EmptyPyGetSetDef_7,[7]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_8 EmptyPyMethodDef_7,[7]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_8 EmptyPyMemberDef_7,[7]={NULL,0,0,0,NULL}
#define DESIG_9(v) DESIG_8(v),[8]=v
#define DESIG__9(n,v) DESIG__8(n,v),[8]={DESIG_##n(v)}
#define EmptyPyGetSetDef_9 EmptyPyGetSetDef_8,[8]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_9 EmptyPyMethodDef_8,[8]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_9 EmptyPyMemberDef_8,[8]={NULL,0,0,0,NULL}
#define DESIG_10(v) DESIG_9(v),[9]=v
#define DESIG__10(n,v) DESIG__9(n,v),[9]={DESIG_##n(v)}
#define EmptyPyGetSetDef_10 EmptyPyGetSetDef_9,[9]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_10 EmptyPyMethodDef_9,[9]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_10 EmptyPyMemberDef_9,[9]={NULL,0,0,0,NULL}
#define DESIG_11(v) DESIG_10(v),[10]=v
#define DESIG__11(n,v) DESIG__10(n,v),[10]={DESIG_##n(v)}
#define EmptyPyGetSetDef_11 EmptyPyGetSetDef_10,[10]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_11 EmptyPyMethodDef_10,[10]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_11 EmptyPyMemberDef_10,[10]={NULL,0,0,0,NULL}
#define DESIG_12(v) DESIG_11(v),[11]=v
#define DESIG__12(n,v) DESIG__11(n,v),[11]={DESIG_##n(v)}
#define EmptyPyGetSetDef_12 EmptyPyGetSetDef_11,[11]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_12 EmptyPyMethodDef_11,[11]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_12 EmptyPyMemberDef_11,[11]={NULL,0,0,0,NULL}
#define DESIG_13(v) DESIG_12(v),[12]=v
#define DESIG__13(n,v) DESIG__12(n,v),[12]={DESIG_##n(v)}
#define EmptyPyGetSetDef_13 EmptyPyGetSetDef_12,[12]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_13 EmptyPyMethodDef_12,[12]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_13 EmptyPyMemberDef_12,[12]={NULL,0,0,0,NULL}
#define DESIG_14(v) DESIG_13(v),[13]=v
#define DESIG__14(n,v) DESIG__13(n,v),[13]={DESIG_##n(v)}
#define EmptyPyGetSetDef_14 EmptyPyGetSetDef_13,[13]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_14 EmptyPyMethodDef_13,[13]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_14 EmptyPyMemberDef_13,[13]={NULL,0,0,0,NULL}
#define DESIG_15(v) DESIG_14(v),[14]=v
#define DESIG__15(n,v) DESIG__14(n,v),[14]={DESIG_##n(v)}
#define EmptyPyGetSetDef_15 EmptyPyGetSetDef_14,[14]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_15 EmptyPyMethodDef_14,[14]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_15 EmptyPyMemberDef_14,[14]={NULL,0,0,0,NULL}
#define DESIG_16(v) DESIG_15(v),[15]=v
#define DESIG__16(n,v) DESIG__15(n,v),[15]={DESIG_##n(v)}
#define EmptyPyGetSetDef_16 EmptyPyGetSetDef_15,[15]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_16 EmptyPyMethodDef_15,[15]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_16 EmptyPyMemberDef_15,[15]={NULL,0,0,0,NULL}
#define DESIG_17(v) DESIG_16(v),[16]=v
#define DESIG__17(n,v) DESIG__16(n,v),[16]={DESIG_##n(v)}
#define EmptyPyGetSetDef_17 EmptyPyGetSetDef_16,[16]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_17 EmptyPyMethodDef_16,[16]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_17 EmptyPyMemberDef_16,[16]={NULL,0,0,0,NULL}
#define DESIG_18(v) DESIG_17(v),[17]=v
#define DESIG__18(n,v) DESIG__17(n,v),[17]={DESIG_##n(v)}
#define EmptyPyGetSetDef_18 EmptyPyGetSetDef_17,[17]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_18 EmptyPyMethodDef_17,[17]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_18 EmptyPyMemberDef_17,[17]={NULL,0,0,0,NULL}
#define DESIG_19(v) DESIG_18(v),[18]=v
#define DESIG__19(n,v) DESIG__18(n,v),[18]={DESIG_##n(v)}
#define EmptyPyGetSetDef_19 EmptyPyGetSetDef_18,[18]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_19 EmptyPyMethodDef_18,[18]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_19 EmptyPyMemberDef_18,[18]={NULL,0,0,0,NULL}
#define DESIG_20(v) DESIG_19(v),[19]=v
#define DESIG__20(n,v) DESIG__19(n,v),[19]={DESIG_##n(v)}
#define EmptyPyGetSetDef_20 EmptyPyGetSetDef_19,[19]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_20 EmptyPyMethodDef_19,[19]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_20 EmptyPyMemberDef_19,[19]={NULL,0,0,0,NULL}
#define DESIG_21(v) DESIG_20(v),[20]=v
#define DESIG__21(n,v) DESIG__20(n,v),[20]={DESIG_##n(v)}
#define EmptyPyGetSetDef_21 EmptyPyGetSetDef_20,[20]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_21 EmptyPyMethodDef_20,[20]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_21 EmptyPyMemberDef_20,[20]={NULL,0,0,0,NULL}
#define DESIG_22(v) DESIG_21(v),[21]=v
#define DESIG__22(n,v) DESIG__21(n,v),[21]={DESIG_##n(v)}
#define EmptyPyGetSetDef_22 EmptyPyGetSetDef_21,[21]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_22 EmptyPyMethodDef_21,[21]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_22 EmptyPyMemberDef_21,[21]={NULL,0,0,0,NULL}
#define DESIG_23(v) DESIG_22(v),[22]=v
#define DESIG__23(n,v) DESIG__22(n,v),[22]={DESIG_##n(v)}
#define EmptyPyGetSetDef_23 EmptyPyGetSetDef_22,[22]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_23 EmptyPyMethodDef_22,[22]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_23 EmptyPyMemberDef_22,[22]={NULL,0,0,0,NULL}
#define DESIG_24(v) DESIG_23(v),[23]=v
#define DESIG__24(n,v) DESIG__23(n,v),[23]={DESIG_##n(v)}
#define EmptyPyGetSetDef_24 EmptyPyGetSetDef_23,[23]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_24 EmptyPyMethodDef_23,[23]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_24 EmptyPyMemberDef_23,[23]={NULL,0,0,0,NULL}
#define DESIG_25(v) DESIG_24(v),[24]=v
#define DESIG__25(n,v) DESIG__24(n,v),[24]={DESIG_##n(v)}
#define EmptyPyGetSetDef_25 EmptyPyGetSetDef_24,[24]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_25 EmptyPyMethodDef_24,[24]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_25 EmptyPyMemberDef_24,[24]={NULL,0,0,0,NULL}
#define DESIG_26(v) DESIG_25(v),[25]=v
#define DESIG__26(n,v) DESIG__25(n,v),[25]={DESIG_##n(v)}
#define EmptyPyGetSetDef_26 EmptyPyGetSetDef_25,[25]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_26 EmptyPyMethodDef_25,[25]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_26 EmptyPyMemberDef_25,[25]={NULL,0,0,0,NULL}
#define DESIG_27(v) DESIG_26(v),[26]=v
#define DESIG__27(n,v) DESIG__26(n,v),[26]={DESIG_##n(v)}
#define EmptyPyGetSetDef_27 EmptyPyGetSetDef_26,[26]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_27 EmptyPyMethodDef_26,[26]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_27 EmptyPyMemberDef_26,[26]={NULL,0,0,0,NULL}
#define DESIG_28(v) DESIG_27(v),[27]=v
#define DESIG__28(n,v) DESIG__27(n,v),[27]={DESIG_##n(v)}
#define EmptyPyGetSetDef_28 EmptyPyGetSetDef_27,[27]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_28 EmptyPyMethodDef_27,[27]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_28 EmptyPyMemberDef_27,[27]={NULL,0,0,0,NULL}
#define DESIG_29(v) DESIG_28(v),[28]=v
#define DESIG__29(n,v) DESIG__28(n,v),[28]={DESIG_##n(v)}
#define EmptyPyGetSetDef_29 EmptyPyGetSetDef_28,[28]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_29 EmptyPyMethodDef_28,[28]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_29 EmptyPyMemberDef_28,[28]={NULL,0,0,0,NULL}
#define DESIG_30(v) DESIG_29(v),[29]=v
#define DESIG__30(n,v) DESIG__29(n,v),[29]={DESIG_##n(v)}
#define EmptyPyGetSetDef_30 EmptyPyGetSetDef_29,[29]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_30 EmptyPyMethodDef_29,[29]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_30 EmptyPyMemberDef_29,[29]={NULL,0,0,0,NULL}
#define DESIG_31(v) DESIG_30(v),[30]=v
#define DESIG__31(n,v) DESIG__30(n,v),[30]={DESIG_##n(v)}
#define EmptyPyGetSetDef_31 EmptyPyGetSetDef_30,[30]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_31 EmptyPyMethodDef_30,[30]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_31 EmptyPyMemberDef_30,[30]={NULL,0,0,0,NULL}
#define DESIG_32(v) DESIG_31(v),[31]=v
#define DESIG__32(n,v) DESIG__31(n,v),[31]={DESIG_##n(v)}
#define EmptyPyGetSetDef_32 EmptyPyGetSetDef_31,[31]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_32 EmptyPyMethodDef_31,[31]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_32 EmptyPyMemberDef_31,[31]={NULL,0,0,0,NULL}
#define DESIG_33(v) DESIG_32(v),[32]=v
#define DESIG__33(n,v) DESIG__32(n,v),[32]={DESIG_##n(v)}
#define EmptyPyGetSetDef_33 EmptyPyGetSetDef_32,[32]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_33 EmptyPyMethodDef_32,[32]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_33 EmptyPyMemberDef_32,[32]={NULL,0,0,0,NULL}
#define DESIG_34(v) DESIG_33(v),[33]=v
#define DESIG__34(n,v) DESIG__33(n,v),[33]={DESIG_##n(v)}
#define EmptyPyGetSetDef_34 EmptyPyGetSetDef_33,[33]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_34 EmptyPyMethodDef_33,[33]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_34 EmptyPyMemberDef_33,[33]={NULL,0,0,0,NULL}
#define DESIG_35(v) DESIG_34(v),[34]=v
#define DESIG__35(n,v) DESIG__34(n,v),[34]={DESIG_##n(v)}
#define EmptyPyGetSetDef_35 EmptyPyGetSetDef_34,[34]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_35 EmptyPyMethodDef_34,[34]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_35 EmptyPyMemberDef_34,[34]={NULL,0,0,0,NULL}
#define DESIG_36(v) DESIG_35(v),[35]=v
#define DESIG__36(n,v) DESIG__35(n,v),[35]={DESIG_##n(v)}
#define EmptyPyGetSetDef_36 EmptyPyGetSetDef_35,[35]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_36 EmptyPyMethodDef_35,[35]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_36 EmptyPyMemberDef_35,[35]={NULL,0,0,0,NULL}
#define DESIG_37(v) DESIG_36(v),[36]=v
#define DESIG__37(n,v) DESIG__36(n,v),[36]={DESIG_##n(v)}
#define EmptyPyGetSetDef_37 EmptyPyGetSetDef_36,[36]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_37 EmptyPyMethodDef_36,[36]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_37 EmptyPyMemberDef_36,[36]={NULL,0,0,0,NULL}
#define DESIG_38(v) DESIG_37(v),[37]=v
#define DESIG__38(n,v) DESIG__37(n,v),[37]={DESIG_##n(v)}
#define EmptyPyGetSetDef_38 EmptyPyGetSetDef_37,[37]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_38 EmptyPyMethodDef_37,[37]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_38 EmptyPyMemberDef_37,[37]={NULL,0,0,0,NULL}
#define DESIG_39(v) DESIG_38(v),[38]=v
#define DESIG__39(n,v) DESIG__38(n,v),[38]={DESIG_##n(v)}
#define EmptyPyGetSetDef_39 EmptyPyGetSetDef_38,[38]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_39 EmptyPyMethodDef_38,[38]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_39 EmptyPyMemberDef_38,[38]={NULL,0,0,0,NULL}
#define DESIG_40(v) DESIG_39(v),[39]=v
#define DESIG__40(n,v) DESIG__39(n,v),[39]={DESIG_##n(v)}
#define EmptyPyGetSetDef_40 EmptyPyGetSetDef_39,[39]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_40 EmptyPyMethodDef_39,[39]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_40 EmptyPyMemberDef_39,[39]={NULL,0,0,0,NULL}
#define DESIG_41(v) DESIG_40(v),[40]=v
#define DESIG__41(n,v) DESIG__40(n,v),[40]={DESIG_##n(v)}
#define EmptyPyGetSetDef_41 EmptyPyGetSetDef_40,[40]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_41 EmptyPyMethodDef_40,[40]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_41 EmptyPyMemberDef_40,[40]={NULL,0,0,0,NULL}
#define DESIG_42(v) DESIG_41(v),[41]=v
#define DESIG__42(n,v) DESIG__41(n,v),[41]={DESIG_##n(v)}
#define EmptyPyGetSetDef_42 EmptyPyGetSetDef_41,[41]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_42 EmptyPyMethodDef_41,[41]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_42 EmptyPyMemberDef_41,[41]={NULL,0,0,0,NULL}
#define DESIG_43(v) DESIG_42(v),[42]=v
#define DESIG__43(n,v) DESIG__42(n,v),[42]={DESIG_##n(v)}
#define EmptyPyGetSetDef_43 EmptyPyGetSetDef_42,[42]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_43 EmptyPyMethodDef_42,[42]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_43 EmptyPyMemberDef_42,[42]={NULL,0,0,0,NULL}
#define DESIG_44(v) DESIG_43(v),[43]=v
#define DESIG__44(n,v) DESIG__43(n,v),[43]={DESIG_##n(v)}
#define EmptyPyGetSetDef_44 EmptyPyGetSetDef_43,[43]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_44 EmptyPyMethodDef_43,[43]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_44 EmptyPyMemberDef_43,[43]={NULL,0,0,0,NULL}
#define DESIG_45(v) DESIG_44(v),[44]=v
#define DESIG__45(n,v) DESIG__44(n,v),[44]={DESIG_##n(v)}
#define EmptyPyGetSetDef_45 EmptyPyGetSetDef_44,[44]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_45 EmptyPyMethodDef_44,[44]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_45 EmptyPyMemberDef_44,[44]={NULL,0,0,0,NULL}
#define DESIG_46(v) DESIG_45(v),[45]=v
#define DESIG__46(n,v) DESIG__45(n,v),[45]={DESIG_##n(v)}
#define EmptyPyGetSetDef_46 EmptyPyGetSetDef_45,[45]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_46 EmptyPyMethodDef_45,[45]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_46 EmptyPyMemberDef_45,[45]={NULL,0,0,0,NULL}
#define DESIG_47(v) DESIG_46(v),[46]=v
#define DESIG__47(n,v) DESIG__46(n,v),[46]={DESIG_##n(v)}
#define EmptyPyGetSetDef_47 EmptyPyGetSetDef_46,[46]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_47 EmptyPyMethodDef_46,[46]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_47 EmptyPyMemberDef_46,[46]={NULL,0,0,0,NULL}
#define DESIG_48(v) DESIG_47(v),[47]=v
#define DESIG__48(n,v) DESIG__47(n,v),[47]={DESIG_##n(v)}
#define EmptyPyGetSetDef_48 EmptyPyGetSetDef_47,[47]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_48 EmptyPyMethodDef_47,[47]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_48 EmptyPyMemberDef_47,[47]={NULL,0,0,0,NULL}
#define DESIG_49(v) DESIG_48(v),[48]=v
#define DESIG__49(n,v) DESIG__48(n,v),[48]={DESIG_##n(v)}
#define EmptyPyGetSetDef_49 EmptyPyGetSetDef_48,[48]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_49 EmptyPyMethodDef_48,[48]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_49 EmptyPyMemberDef_48,[48]={NULL,0,0,0,NULL}
#define DESIG_50(v) DESIG_49(v),[49]=v
#define DESIG__50(n,v) DESIG__49(n,v),[49]={DESIG_##n(v)}
#define EmptyPyGetSetDef_50 EmptyPyGetSetDef_49,[49]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_50 EmptyPyMethodDef_49,[49]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_50 EmptyPyMemberDef_49,[49]={NULL,0,0,0,NULL}
#define DESIG_51(v) DESIG_50(v),[50]=v
#define DESIG__51(n,v) DESIG__50(n,v),[50]={DESIG_##n(v)}
#define EmptyPyGetSetDef_51 EmptyPyGetSetDef_50,[50]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_51 EmptyPyMethodDef_50,[50]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_51 EmptyPyMemberDef_50,[50]={NULL,0,0,0,NULL}
#define DESIG_52(v) DESIG_51(v),[51]=v
#define DESIG__52(n,v) DESIG__51(n,v),[51]={DESIG_##n(v)}
#define EmptyPyGetSetDef_52 EmptyPyGetSetDef_51,[51]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_52 EmptyPyMethodDef_51,[51]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_52 EmptyPyMemberDef_51,[51]={NULL,0,0,0,NULL}
#define DESIG_53(v) DESIG_52(v),[52]=v
#define DESIG__53(n,v) DESIG__52(n,v),[52]={DESIG_##n(v)}
#define EmptyPyGetSetDef_53 EmptyPyGetSetDef_52,[52]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_53 EmptyPyMethodDef_52,[52]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_53 EmptyPyMemberDef_52,[52]={NULL,0,0,0,NULL}
#define DESIG_54(v) DESIG_53(v),[53]=v
#define DESIG__54(n,v) DESIG__53(n,v),[53]={DESIG_##n(v)}
#define EmptyPyGetSetDef_54 EmptyPyGetSetDef_53,[53]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_54 EmptyPyMethodDef_53,[53]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_54 EmptyPyMemberDef_53,[53]={NULL,0,0,0,NULL}
#define DESIG_55(v) DESIG_54(v),[54]=v
#define DESIG__55(n,v) DESIG__54(n,v),[54]={DESIG_##n(v)}
#define EmptyPyGetSetDef_55 EmptyPyGetSetDef_54,[54]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_55 EmptyPyMethodDef_54,[54]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_55 EmptyPyMemberDef_54,[54]={NULL,0,0,0,NULL}
#define DESIG_56(v) DESIG_55(v),[55]=v
#define DESIG__56(n,v) DESIG__55(n,v),[55]={DESIG_##n(v)}
#define EmptyPyGetSetDef_56 EmptyPyGetSetDef_55,[55]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_56 EmptyPyMethodDef_55,[55]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_56 EmptyPyMemberDef_55,[55]={NULL,0,0,0,NULL}
#define DESIG_57(v) DESIG_56(v),[56]=v
#define DESIG__57(n,v) DESIG__56(n,v),[56]={DESIG_##n(v)}
#define EmptyPyGetSetDef_57 EmptyPyGetSetDef_56,[56]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_57 EmptyPyMethodDef_56,[56]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_57 EmptyPyMemberDef_56,[56]={NULL,0,0,0,NULL}
#define DESIG_58(v) DESIG_57(v),[57]=v
#define DESIG__58(n,v) DESIG__57(n,v),[57]={DESIG_##n(v)}
#define EmptyPyGetSetDef_58 EmptyPyGetSetDef_57,[57]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_58 EmptyPyMethodDef_57,[57]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_58 EmptyPyMemberDef_57,[57]={NULL,0,0,0,NULL}
#define DESIG_59(v) DESIG_58(v),[58]=v
#define DESIG__59(n,v) DESIG__58(n,v),[58]={DESIG_##n(v)}
#define EmptyPyGetSetDef_59 EmptyPyGetSetDef_58,[58]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_59 EmptyPyMethodDef_58,[58]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_59 EmptyPyMemberDef_58,[58]={NULL,0,0,0,NULL}
#define DESIG_60(v) DESIG_59(v),[59]=v
#define DESIG__60(n,v) DESIG__59(n,v),[59]={DESIG_##n(v)}
#define EmptyPyGetSetDef_60 EmptyPyGetSetDef_59,[59]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_60 EmptyPyMethodDef_59,[59]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_60 EmptyPyMemberDef_59,[59]={NULL,0,0,0,NULL}
#define DESIG_61(v) DESIG_60(v),[60]=v
#define DESIG__61(n,v) DESIG__60(n,v),[60]={DESIG_##n(v)}
#define EmptyPyGetSetDef_61 EmptyPyGetSetDef_60,[60]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_61 EmptyPyMethodDef_60,[60]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_61 EmptyPyMemberDef_60,[60]={NULL,0,0,0,NULL}
#define DESIG_62(v) DESIG_61(v),[61]=v
#define DESIG__62(n,v) DESIG__61(n,v),[61]={DESIG_##n(v)}
#define EmptyPyGetSetDef_62 EmptyPyGetSetDef_61,[61]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_62 EmptyPyMethodDef_61,[61]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_62 EmptyPyMemberDef_61,[61]={NULL,0,0,0,NULL}
#define DESIG_63(v) DESIG_62(v),[62]=v
#define DESIG__63(n,v) DESIG__62(n,v),[62]={DESIG_##n(v)}
#define EmptyPyGetSetDef_63 EmptyPyGetSetDef_62,[62]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_63 EmptyPyMethodDef_62,[62]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_63 EmptyPyMemberDef_62,[62]={NULL,0,0,0,NULL}
#define DESIG_64(v) DESIG_63(v),[63]=v
#define DESIG__64(n,v) DESIG__63(n,v),[63]={DESIG_##n(v)}
#define EmptyPyGetSetDef_64 EmptyPyGetSetDef_63,[63]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_64 EmptyPyMethodDef_63,[63]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_64 EmptyPyMemberDef_63,[63]={NULL,0,0,0,NULL}
#define DESIG_65(v) DESIG_64(v),[64]=v
#define DESIG__65(n,v) DESIG__64(n,v),[64]={DESIG_##n(v)}
#define EmptyPyGetSetDef_65 EmptyPyGetSetDef_64,[64]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_65 EmptyPyMethodDef_64,[64]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_65 EmptyPyMemberDef_64,[64]={NULL,0,0,0,NULL}
#define DESIG_66(v) DESIG_65(v),[65]=v
#define DESIG__66(n,v) DESIG__65(n,v),[65]={DESIG_##n(v)}
#define EmptyPyGetSetDef_66 EmptyPyGetSetDef_65,[65]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_66 EmptyPyMethodDef_65,[65]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_66 EmptyPyMemberDef_65,[65]={NULL,0,0,0,NULL}
#define DESIG_67(v) DESIG_66(v),[66]=v
#define DESIG__67(n,v) DESIG__66(n,v),[66]={DESIG_##n(v)}
#define EmptyPyGetSetDef_67 EmptyPyGetSetDef_66,[66]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_67 EmptyPyMethodDef_66,[66]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_67 EmptyPyMemberDef_66,[66]={NULL,0,0,0,NULL}
#define DESIG_68(v) DESIG_67(v),[67]=v
#define DESIG__68(n,v) DESIG__67(n,v),[67]={DESIG_##n(v)}
#define EmptyPyGetSetDef_68 EmptyPyGetSetDef_67,[67]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_68 EmptyPyMethodDef_67,[67]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_68 EmptyPyMemberDef_67,[67]={NULL,0,0,0,NULL}
#define DESIG_69(v) DESIG_68(v),[68]=v
#define DESIG__69(n,v) DESIG__68(n,v),[68]={DESIG_##n(v)}
#define EmptyPyGetSetDef_69 EmptyPyGetSetDef_68,[68]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_69 EmptyPyMethodDef_68,[68]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_69 EmptyPyMemberDef_68,[68]={NULL,0,0,0,NULL}
#define DESIG_70(v) DESIG_69(v),[69]=v
#define DESIG__70(n,v) DESIG__69(n,v),[69]={DESIG_##n(v)}
#define EmptyPyGetSetDef_70 EmptyPyGetSetDef_69,[69]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_70 EmptyPyMethodDef_69,[69]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_70 EmptyPyMemberDef_69,[69]={NULL,0,0,0,NULL}
#define DESIG_71(v) DESIG_70(v),[70]=v
#define DESIG__71(n,v) DESIG__70(n,v),[70]={DESIG_##n(v)}
#define EmptyPyGetSetDef_71 EmptyPyGetSetDef_70,[70]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_71 EmptyPyMethodDef_70,[70]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_71 EmptyPyMemberDef_70,[70]={NULL,0,0,0,NULL}
#define DESIG_72(v) DESIG_71(v),[71]=v
#define DESIG__72(n,v) DESIG__71(n,v),[71]={DESIG_##n(v)}
#define EmptyPyGetSetDef_72 EmptyPyGetSetDef_71,[71]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_72 EmptyPyMethodDef_71,[71]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_72 EmptyPyMemberDef_71,[71]={NULL,0,0,0,NULL}
#define DESIG_73(v) DESIG_72(v),[72]=v
#define DESIG__73(n,v) DESIG__72(n,v),[72]={DESIG_##n(v)}
#define EmptyPyGetSetDef_73 EmptyPyGetSetDef_72,[72]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_73 EmptyPyMethodDef_72,[72]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_73 EmptyPyMemberDef_72,[72]={NULL,0,0,0,NULL}
#define DESIG_74(v) DESIG_73(v),[73]=v
#define DESIG__74(n,v) DESIG__73(n,v),[73]={DESIG_##n(v)}
#define EmptyPyGetSetDef_74 EmptyPyGetSetDef_73,[73]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_74 EmptyPyMethodDef_73,[73]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_74 EmptyPyMemberDef_73,[73]={NULL,0,0,0,NULL}
#define DESIG_75(v) DESIG_74(v),[74]=v
#define DESIG__75(n,v) DESIG__74(n,v),[74]={DESIG_##n(v)}
#define EmptyPyGetSetDef_75 EmptyPyGetSetDef_74,[74]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_75 EmptyPyMethodDef_74,[74]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_75 EmptyPyMemberDef_74,[74]={NULL,0,0,0,NULL}
#define DESIG_76(v) DESIG_75(v),[75]=v
#define DESIG__76(n,v) DESIG__75(n,v),[75]={DESIG_##n(v)}
#define EmptyPyGetSetDef_76 EmptyPyGetSetDef_75,[75]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_76 EmptyPyMethodDef_75,[75]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_76 EmptyPyMemberDef_75,[75]={NULL,0,0,0,NULL}
#define DESIG_77(v) DESIG_76(v),[76]=v
#define DESIG__77(n,v) DESIG__76(n,v),[76]={DESIG_##n(v)}
#define EmptyPyGetSetDef_77 EmptyPyGetSetDef_76,[76]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_77 EmptyPyMethodDef_76,[76]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_77 EmptyPyMemberDef_76,[76]={NULL,0,0,0,NULL}
#define DESIG_78(v) DESIG_77(v),[77]=v
#define DESIG__78(n,v) DESIG__77(n,v),[77]={DESIG_##n(v)}
#define EmptyPyGetSetDef_78 EmptyPyGetSetDef_77,[77]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_78 EmptyPyMethodDef_77,[77]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_78 EmptyPyMemberDef_77,[77]={NULL,0,0,0,NULL}
#define DESIG_79(v) DESIG_78(v),[78]=v
#define DESIG__79(n,v) DESIG__78(n,v),[78]={DESIG_##n(v)}
#define EmptyPyGetSetDef_79 EmptyPyGetSetDef_78,[78]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_79 EmptyPyMethodDef_78,[78]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_79 EmptyPyMemberDef_78,[78]={NULL,0,0,0,NULL}
#define DESIG_80(v) DESIG_79(v),[79]=v
#define DESIG__80(n,v) DESIG__79(n,v),[79]={DESIG_##n(v)}
#define EmptyPyGetSetDef_80 EmptyPyGetSetDef_79,[79]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_80 EmptyPyMethodDef_79,[79]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_80 EmptyPyMemberDef_79,[79]={NULL,0,0,0,NULL}
#define DESIG_81(v) DESIG_80(v),[80]=v
#define DESIG__81(n,v) DESIG__80(n,v),[80]={DESIG_##n(v)}
#define EmptyPyGetSetDef_81 EmptyPyGetSetDef_80,[80]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_81 EmptyPyMethodDef_80,[80]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_81 EmptyPyMemberDef_80,[80]={NULL,0,0,0,NULL}
#define DESIG_82(v) DESIG_81(v),[81]=v
#define DESIG__82(n,v) DESIG__81(n,v),[81]={DESIG_##n(v)}
#define EmptyPyGetSetDef_82 EmptyPyGetSetDef_81,[81]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_82 EmptyPyMethodDef_81,[81]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_82 EmptyPyMemberDef_81,[81]={NULL,0,0,0,NULL}
#define DESIG_83(v) DESIG_82(v),[82]=v
#define DESIG__83(n,v) DESIG__82(n,v),[82]={DESIG_##n(v)}
#define EmptyPyGetSetDef_83 EmptyPyGetSetDef_82,[82]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_83 EmptyPyMethodDef_82,[82]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_83 EmptyPyMemberDef_82,[82]={NULL,0,0,0,NULL}
#define DESIG_84(v) DESIG_83(v),[83]=v
#define DESIG__84(n,v) DESIG__83(n,v),[83]={DESIG_##n(v)}
#define EmptyPyGetSetDef_84 EmptyPyGetSetDef_83,[83]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_84 EmptyPyMethodDef_83,[83]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_84 EmptyPyMemberDef_83,[83]={NULL,0,0,0,NULL}
#define DESIG_85(v) DESIG_84(v),[84]=v
#define DESIG__85(n,v) DESIG__84(n,v),[84]={DESIG_##n(v)}
#define EmptyPyGetSetDef_85 EmptyPyGetSetDef_84,[84]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_85 EmptyPyMethodDef_84,[84]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_85 EmptyPyMemberDef_84,[84]={NULL,0,0,0,NULL}
#define DESIG_86(v) DESIG_85(v),[85]=v
#define DESIG__86(n,v) DESIG__85(n,v),[85]={DESIG_##n(v)}
#define EmptyPyGetSetDef_86 EmptyPyGetSetDef_85,[85]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_86 EmptyPyMethodDef_85,[85]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_86 EmptyPyMemberDef_85,[85]={NULL,0,0,0,NULL}
#define DESIG_87(v) DESIG_86(v),[86]=v
#define DESIG__87(n,v) DESIG__86(n,v),[86]={DESIG_##n(v)}
#define EmptyPyGetSetDef_87 EmptyPyGetSetDef_86,[86]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_87 EmptyPyMethodDef_86,[86]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_87 EmptyPyMemberDef_86,[86]={NULL,0,0,0,NULL}
#define DESIG_88(v) DESIG_87(v),[87]=v
#define DESIG__88(n,v) DESIG__87(n,v),[87]={DESIG_##n(v)}
#define EmptyPyGetSetDef_88 EmptyPyGetSetDef_87,[87]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_88 EmptyPyMethodDef_87,[87]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_88 EmptyPyMemberDef_87,[87]={NULL,0,0,0,NULL}
#define DESIG_89(v) DESIG_88(v),[88]=v
#define DESIG__89(n,v) DESIG__88(n,v),[88]={DESIG_##n(v)}
#define EmptyPyGetSetDef_89 EmptyPyGetSetDef_88,[88]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_89 EmptyPyMethodDef_88,[88]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_89 EmptyPyMemberDef_88,[88]={NULL,0,0,0,NULL}
#define DESIG_90(v) DESIG_89(v),[89]=v
#define DESIG__90(n,v) DESIG__89(n,v),[89]={DESIG_##n(v)}
#define EmptyPyGetSetDef_90 EmptyPyGetSetDef_89,[89]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_90 EmptyPyMethodDef_89,[89]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_90 EmptyPyMemberDef_89,[89]={NULL,0,0,0,NULL}
#define DESIG_91(v) DESIG_90(v),[90]=v
#define DESIG__91(n,v) DESIG__90(n,v),[90]={DESIG_##n(v)}
#define EmptyPyGetSetDef_91 EmptyPyGetSetDef_90,[90]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_91 EmptyPyMethodDef_90,[90]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_91 EmptyPyMemberDef_90,[90]={NULL,0,0,0,NULL}
#define DESIG_92(v) DESIG_91(v),[91]=v
#define DESIG__92(n,v) DESIG__91(n,v),[91]={DESIG_##n(v)}
#define EmptyPyGetSetDef_92 EmptyPyGetSetDef_91,[91]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_92 EmptyPyMethodDef_91,[91]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_92 EmptyPyMemberDef_91,[91]={NULL,0,0,0,NULL}
#define DESIG_93(v) DESIG_92(v),[92]=v
#define DESIG__93(n,v) DESIG__92(n,v),[92]={DESIG_##n(v)}
#define EmptyPyGetSetDef_93 EmptyPyGetSetDef_92,[92]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_93 EmptyPyMethodDef_92,[92]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_93 EmptyPyMemberDef_92,[92]={NULL,0,0,0,NULL}
#define DESIG_94(v) DESIG_93(v),[93]=v
#define DESIG__94(n,v) DESIG__93(n,v),[93]={DESIG_##n(v)}
#define EmptyPyGetSetDef_94 EmptyPyGetSetDef_93,[93]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_94 EmptyPyMethodDef_93,[93]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_94 EmptyPyMemberDef_93,[93]={NULL,0,0,0,NULL}
#define DESIG_95(v) DESIG_94(v),[94]=v
#define DESIG__95(n,v) DESIG__94(n,v),[94]={DESIG_##n(v)}
#define EmptyPyGetSetDef_95 EmptyPyGetSetDef_94,[94]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_95 EmptyPyMethodDef_94,[94]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_95 EmptyPyMemberDef_94,[94]={NULL,0,0,0,NULL}
#define DESIG_96(v) DESIG_95(v),[95]=v
#define DESIG__96(n,v) DESIG__95(n,v),[95]={DESIG_##n(v)}
#define EmptyPyGetSetDef_96 EmptyPyGetSetDef_95,[95]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_96 EmptyPyMethodDef_95,[95]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_96 EmptyPyMemberDef_95,[95]={NULL,0,0,0,NULL}
#define DESIG_97(v) DESIG_96(v),[96]=v
#define DESIG__97(n,v) DESIG__96(n,v),[96]={DESIG_##n(v)}
#define EmptyPyGetSetDef_97 EmptyPyGetSetDef_96,[96]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_97 EmptyPyMethodDef_96,[96]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_97 EmptyPyMemberDef_96,[96]={NULL,0,0,0,NULL}
#define DESIG_98(v) DESIG_97(v),[97]=v
#define DESIG__98(n,v) DESIG__97(n,v),[97]={DESIG_##n(v)}
#define EmptyPyGetSetDef_98 EmptyPyGetSetDef_97,[97]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_98 EmptyPyMethodDef_97,[97]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_98 EmptyPyMemberDef_97,[97]={NULL,0,0,0,NULL}
#define DESIG_99(v) DESIG_98(v),[98]=v
#define DESIG__99(n,v) DESIG__98(n,v),[98]={DESIG_##n(v)}
#define EmptyPyGetSetDef_99 EmptyPyGetSetDef_98,[98]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_99 EmptyPyMethodDef_98,[98]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_99 EmptyPyMemberDef_98,[98]={NULL,0,0,0,NULL}
#define DESIG_100(v) DESIG_99(v),[99]=v
#define DESIG__100(n,v) DESIG__99(n,v),[99]={DESIG_##n(v)}
#define EmptyPyGetSetDef_100 EmptyPyGetSetDef_99,[99]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_100 EmptyPyMethodDef_99,[99]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_100 EmptyPyMemberDef_99,[99]={NULL,0,0,0,NULL}
#define DESIG_101(v) DESIG_100(v),[100]=v
#define DESIG__101(n,v) DESIG__100(n,v),[100]={DESIG_##n(v)}
#define EmptyPyGetSetDef_101 EmptyPyGetSetDef_100,[100]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_101 EmptyPyMethodDef_100,[100]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_101 EmptyPyMemberDef_100,[100]={NULL,0,0,0,NULL}
#define DESIG_102(v) DESIG_101(v),[101]=v
#define DESIG__102(n,v) DESIG__101(n,v),[101]={DESIG_##n(v)}
#define EmptyPyGetSetDef_102 EmptyPyGetSetDef_101,[101]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_102 EmptyPyMethodDef_101,[101]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_102 EmptyPyMemberDef_101,[101]={NULL,0,0,0,NULL}
#define DESIG_103(v) DESIG_102(v),[102]=v
#define DESIG__103(n,v) DESIG__102(n,v),[102]={DESIG_##n(v)}
#define EmptyPyGetSetDef_103 EmptyPyGetSetDef_102,[102]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_103 EmptyPyMethodDef_102,[102]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_103 EmptyPyMemberDef_102,[102]={NULL,0,0,0,NULL}
#define DESIG_104(v) DESIG_103(v),[103]=v
#define DESIG__104(n,v) DESIG__103(n,v),[103]={DESIG_##n(v)}
#define EmptyPyGetSetDef_104 EmptyPyGetSetDef_103,[103]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_104 EmptyPyMethodDef_103,[103]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_104 EmptyPyMemberDef_103,[103]={NULL,0,0,0,NULL}
#define DESIG_105(v) DESIG_104(v),[104]=v
#define DESIG__105(n,v) DESIG__104(n,v),[104]={DESIG_##n(v)}
#define EmptyPyGetSetDef_105 EmptyPyGetSetDef_104,[104]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_105 EmptyPyMethodDef_104,[104]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_105 EmptyPyMemberDef_104,[104]={NULL,0,0,0,NULL}
#define DESIG_106(v) DESIG_105(v),[105]=v
#define DESIG__106(n,v) DESIG__105(n,v),[105]={DESIG_##n(v)}
#define EmptyPyGetSetDef_106 EmptyPyGetSetDef_105,[105]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_106 EmptyPyMethodDef_105,[105]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_106 EmptyPyMemberDef_105,[105]={NULL,0,0,0,NULL}
#define DESIG_107(v) DESIG_106(v),[106]=v
#define DESIG__107(n,v) DESIG__106(n,v),[106]={DESIG_##n(v)}
#define EmptyPyGetSetDef_107 EmptyPyGetSetDef_106,[106]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_107 EmptyPyMethodDef_106,[106]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_107 EmptyPyMemberDef_106,[106]={NULL,0,0,0,NULL}
#define DESIG_108(v) DESIG_107(v),[107]=v
#define DESIG__108(n,v) DESIG__107(n,v),[107]={DESIG_##n(v)}
#define EmptyPyGetSetDef_108 EmptyPyGetSetDef_107,[107]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_108 EmptyPyMethodDef_107,[107]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_108 EmptyPyMemberDef_107,[107]={NULL,0,0,0,NULL}
#define DESIG_109(v) DESIG_108(v),[108]=v
#define DESIG__109(n,v) DESIG__108(n,v),[108]={DESIG_##n(v)}
#define EmptyPyGetSetDef_109 EmptyPyGetSetDef_108,[108]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_109 EmptyPyMethodDef_108,[108]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_109 EmptyPyMemberDef_108,[108]={NULL,0,0,0,NULL}
#define DESIG_110(v) DESIG_109(v),[109]=v
#define DESIG__110(n,v) DESIG__109(n,v),[109]={DESIG_##n(v)}
#define EmptyPyGetSetDef_110 EmptyPyGetSetDef_109,[109]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_110 EmptyPyMethodDef_109,[109]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_110 EmptyPyMemberDef_109,[109]={NULL,0,0,0,NULL}
#define DESIG_111(v) DESIG_110(v),[110]=v
#define DESIG__111(n,v) DESIG__110(n,v),[110]={DESIG_##n(v)}
#define EmptyPyGetSetDef_111 EmptyPyGetSetDef_110,[110]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_111 EmptyPyMethodDef_110,[110]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_111 EmptyPyMemberDef_110,[110]={NULL,0,0,0,NULL}
#define DESIG_112(v) DESIG_111(v),[111]=v
#define DESIG__112(n,v) DESIG__111(n,v),[111]={DESIG_##n(v)}
#define EmptyPyGetSetDef_112 EmptyPyGetSetDef_111,[111]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_112 EmptyPyMethodDef_111,[111]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_112 EmptyPyMemberDef_111,[111]={NULL,0,0,0,NULL}
#define DESIG_113(v) DESIG_112(v),[112]=v
#define DESIG__113(n,v) DESIG__112(n,v),[112]={DESIG_##n(v)}
#define EmptyPyGetSetDef_113 EmptyPyGetSetDef_112,[112]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_113 EmptyPyMethodDef_112,[112]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_113 EmptyPyMemberDef_112,[112]={NULL,0,0,0,NULL}
#define DESIG_114(v) DESIG_113(v),[113]=v
#define DESIG__114(n,v) DESIG__113(n,v),[113]={DESIG_##n(v)}
#define EmptyPyGetSetDef_114 EmptyPyGetSetDef_113,[113]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_114 EmptyPyMethodDef_113,[113]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_114 EmptyPyMemberDef_113,[113]={NULL,0,0,0,NULL}
#define DESIG_115(v) DESIG_114(v),[114]=v
#define DESIG__115(n,v) DESIG__114(n,v),[114]={DESIG_##n(v)}
#define EmptyPyGetSetDef_115 EmptyPyGetSetDef_114,[114]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_115 EmptyPyMethodDef_114,[114]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_115 EmptyPyMemberDef_114,[114]={NULL,0,0,0,NULL}
#define DESIG_116(v) DESIG_115(v),[115]=v
#define DESIG__116(n,v) DESIG__115(n,v),[115]={DESIG_##n(v)}
#define EmptyPyGetSetDef_116 EmptyPyGetSetDef_115,[115]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_116 EmptyPyMethodDef_115,[115]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_116 EmptyPyMemberDef_115,[115]={NULL,0,0,0,NULL}
#define DESIG_117(v) DESIG_116(v),[116]=v
#define DESIG__117(n,v) DESIG__116(n,v),[116]={DESIG_##n(v)}
#define EmptyPyGetSetDef_117 EmptyPyGetSetDef_116,[116]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_117 EmptyPyMethodDef_116,[116]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_117 EmptyPyMemberDef_116,[116]={NULL,0,0,0,NULL}
#define DESIG_118(v) DESIG_117(v),[117]=v
#define DESIG__118(n,v) DESIG__117(n,v),[117]={DESIG_##n(v)}
#define EmptyPyGetSetDef_118 EmptyPyGetSetDef_117,[117]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_118 EmptyPyMethodDef_117,[117]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_118 EmptyPyMemberDef_117,[117]={NULL,0,0,0,NULL}
#define DESIG_119(v) DESIG_118(v),[118]=v
#define DESIG__119(n,v) DESIG__118(n,v),[118]={DESIG_##n(v)}
#define EmptyPyGetSetDef_119 EmptyPyGetSetDef_118,[118]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_119 EmptyPyMethodDef_118,[118]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_119 EmptyPyMemberDef_118,[118]={NULL,0,0,0,NULL}
#define DESIG_120(v) DESIG_119(v),[119]=v
#define DESIG__120(n,v) DESIG__119(n,v),[119]={DESIG_##n(v)}
#define EmptyPyGetSetDef_120 EmptyPyGetSetDef_119,[119]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_120 EmptyPyMethodDef_119,[119]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_120 EmptyPyMemberDef_119,[119]={NULL,0,0,0,NULL}
#define DESIG_121(v) DESIG_120(v),[120]=v
#define DESIG__121(n,v) DESIG__120(n,v),[120]={DESIG_##n(v)}
#define EmptyPyGetSetDef_121 EmptyPyGetSetDef_120,[120]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_121 EmptyPyMethodDef_120,[120]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_121 EmptyPyMemberDef_120,[120]={NULL,0,0,0,NULL}
#define DESIG_122(v) DESIG_121(v),[121]=v
#define DESIG__122(n,v) DESIG__121(n,v),[121]={DESIG_##n(v)}
#define EmptyPyGetSetDef_122 EmptyPyGetSetDef_121,[121]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_122 EmptyPyMethodDef_121,[121]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_122 EmptyPyMemberDef_121,[121]={NULL,0,0,0,NULL}
#define DESIG_123(v) DESIG_122(v),[122]=v
#define DESIG__123(n,v) DESIG__122(n,v),[122]={DESIG_##n(v)}
#define EmptyPyGetSetDef_123 EmptyPyGetSetDef_122,[122]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_123 EmptyPyMethodDef_122,[122]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_123 EmptyPyMemberDef_122,[122]={NULL,0,0,0,NULL}
#define DESIG_124(v) DESIG_123(v),[123]=v
#define DESIG__124(n,v) DESIG__123(n,v),[123]={DESIG_##n(v)}
#define EmptyPyGetSetDef_124 EmptyPyGetSetDef_123,[123]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_124 EmptyPyMethodDef_123,[123]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_124 EmptyPyMemberDef_123,[123]={NULL,0,0,0,NULL}
#define DESIG_125(v) DESIG_124(v),[124]=v
#define DESIG__125(n,v) DESIG__124(n,v),[124]={DESIG_##n(v)}
#define EmptyPyGetSetDef_125 EmptyPyGetSetDef_124,[124]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_125 EmptyPyMethodDef_124,[124]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_125 EmptyPyMemberDef_124,[124]={NULL,0,0,0,NULL}
#define DESIG_126(v) DESIG_125(v),[125]=v
#define DESIG__126(n,v) DESIG__125(n,v),[125]={DESIG_##n(v)}
#define EmptyPyGetSetDef_126 EmptyPyGetSetDef_125,[125]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_126 EmptyPyMethodDef_125,[125]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_126 EmptyPyMemberDef_125,[125]={NULL,0,0,0,NULL}
#define DESIG_127(v) DESIG_126(v),[126]=v
#define DESIG__127(n,v) DESIG__126(n,v),[126]={DESIG_##n(v)}
#define EmptyPyGetSetDef_127 EmptyPyGetSetDef_126,[126]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_127 EmptyPyMethodDef_126,[126]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_127 EmptyPyMemberDef_126,[126]={NULL,0,0,0,NULL}
#define DESIG_128(v) DESIG_127(v),[127]=v
#define DESIG__128(n,v) DESIG__127(n,v),[127]={DESIG_##n(v)}
#define EmptyPyGetSetDef_128 EmptyPyGetSetDef_127,[127]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_128 EmptyPyMethodDef_127,[127]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_128 EmptyPyMemberDef_127,[127]={NULL,0,0,0,NULL}
#define DESIG_129(v) DESIG_128(v),[128]=v
#define DESIG__129(n,v) DESIG__128(n,v),[128]={DESIG_##n(v)}
#define EmptyPyGetSetDef_129 EmptyPyGetSetDef_128,[128]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_129 EmptyPyMethodDef_128,[128]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_129 EmptyPyMemberDef_128,[128]={NULL,0,0,0,NULL}
#define DESIG_130(v) DESIG_129(v),[129]=v
#define DESIG__130(n,v) DESIG__129(n,v),[129]={DESIG_##n(v)}
#define EmptyPyGetSetDef_130 EmptyPyGetSetDef_129,[129]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_130 EmptyPyMethodDef_129,[129]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_130 EmptyPyMemberDef_129,[129]={NULL,0,0,0,NULL}
#define DESIG_131(v) DESIG_130(v),[130]=v
#define DESIG__131(n,v) DESIG__130(n,v),[130]={DESIG_##n(v)}
#define EmptyPyGetSetDef_131 EmptyPyGetSetDef_130,[130]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_131 EmptyPyMethodDef_130,[130]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_131 EmptyPyMemberDef_130,[130]={NULL,0,0,0,NULL}
#define DESIG_132(v) DESIG_131(v),[131]=v
#define DESIG__132(n,v) DESIG__131(n,v),[131]={DESIG_##n(v)}
#define EmptyPyGetSetDef_132 EmptyPyGetSetDef_131,[131]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_132 EmptyPyMethodDef_131,[131]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_132 EmptyPyMemberDef_131,[131]={NULL,0,0,0,NULL}
#define DESIG_133(v) DESIG_132(v),[132]=v
#define DESIG__133(n,v) DESIG__132(n,v),[132]={DESIG_##n(v)}
#define EmptyPyGetSetDef_133 EmptyPyGetSetDef_132,[132]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_133 EmptyPyMethodDef_132,[132]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_133 EmptyPyMemberDef_132,[132]={NULL,0,0,0,NULL}
#define DESIG_134(v) DESIG_133(v),[133]=v
#define DESIG__134(n,v) DESIG__133(n,v),[133]={DESIG_##n(v)}
#define EmptyPyGetSetDef_134 EmptyPyGetSetDef_133,[133]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_134 EmptyPyMethodDef_133,[133]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_134 EmptyPyMemberDef_133,[133]={NULL,0,0,0,NULL}
#define DESIG_135(v) DESIG_134(v),[134]=v
#define DESIG__135(n,v) DESIG__134(n,v),[134]={DESIG_##n(v)}
#define EmptyPyGetSetDef_135 EmptyPyGetSetDef_134,[134]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_135 EmptyPyMethodDef_134,[134]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_135 EmptyPyMemberDef_134,[134]={NULL,0,0,0,NULL}
#define DESIG_136(v) DESIG_135(v),[135]=v
#define DESIG__136(n,v) DESIG__135(n,v),[135]={DESIG_##n(v)}
#define EmptyPyGetSetDef_136 EmptyPyGetSetDef_135,[135]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_136 EmptyPyMethodDef_135,[135]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_136 EmptyPyMemberDef_135,[135]={NULL,0,0,0,NULL}
#define DESIG_137(v) DESIG_136(v),[136]=v
#define DESIG__137(n,v) DESIG__136(n,v),[136]={DESIG_##n(v)}
#define EmptyPyGetSetDef_137 EmptyPyGetSetDef_136,[136]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_137 EmptyPyMethodDef_136,[136]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_137 EmptyPyMemberDef_136,[136]={NULL,0,0,0,NULL}
#define DESIG_138(v) DESIG_137(v),[137]=v
#define DESIG__138(n,v) DESIG__137(n,v),[137]={DESIG_##n(v)}
#define EmptyPyGetSetDef_138 EmptyPyGetSetDef_137,[137]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_138 EmptyPyMethodDef_137,[137]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_138 EmptyPyMemberDef_137,[137]={NULL,0,0,0,NULL}
#define DESIG_139(v) DESIG_138(v),[138]=v
#define DESIG__139(n,v) DESIG__138(n,v),[138]={DESIG_##n(v)}
#define EmptyPyGetSetDef_139 EmptyPyGetSetDef_138,[138]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_139 EmptyPyMethodDef_138,[138]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_139 EmptyPyMemberDef_138,[138]={NULL,0,0,0,NULL}
#define DESIG_140(v) DESIG_139(v),[139]=v
#define DESIG__140(n,v) DESIG__139(n,v),[139]={DESIG_##n(v)}
#define EmptyPyGetSetDef_140 EmptyPyGetSetDef_139,[139]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_140 EmptyPyMethodDef_139,[139]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_140 EmptyPyMemberDef_139,[139]={NULL,0,0,0,NULL}
#define DESIG_141(v) DESIG_140(v),[140]=v
#define DESIG__141(n,v) DESIG__140(n,v),[140]={DESIG_##n(v)}
#define EmptyPyGetSetDef_141 EmptyPyGetSetDef_140,[140]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_141 EmptyPyMethodDef_140,[140]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_141 EmptyPyMemberDef_140,[140]={NULL,0,0,0,NULL}
#define DESIG_142(v) DESIG_141(v),[141]=v
#define DESIG__142(n,v) DESIG__141(n,v),[141]={DESIG_##n(v)}
#define EmptyPyGetSetDef_142 EmptyPyGetSetDef_141,[141]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_142 EmptyPyMethodDef_141,[141]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_142 EmptyPyMemberDef_141,[141]={NULL,0,0,0,NULL}
#define DESIG_143(v) DESIG_142(v),[142]=v
#define DESIG__143(n,v) DESIG__142(n,v),[142]={DESIG_##n(v)}
#define EmptyPyGetSetDef_143 EmptyPyGetSetDef_142,[142]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_143 EmptyPyMethodDef_142,[142]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_143 EmptyPyMemberDef_142,[142]={NULL,0,0,0,NULL}
#define DESIG_144(v) DESIG_143(v),[143]=v
#define DESIG__144(n,v) DESIG__143(n,v),[143]={DESIG_##n(v)}
#define EmptyPyGetSetDef_144 EmptyPyGetSetDef_143,[143]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_144 EmptyPyMethodDef_143,[143]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_144 EmptyPyMemberDef_143,[143]={NULL,0,0,0,NULL}
#define DESIG_145(v) DESIG_144(v),[144]=v
#define DESIG__145(n,v) DESIG__144(n,v),[144]={DESIG_##n(v)}
#define EmptyPyGetSetDef_145 EmptyPyGetSetDef_144,[144]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_145 EmptyPyMethodDef_144,[144]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_145 EmptyPyMemberDef_144,[144]={NULL,0,0,0,NULL}
#define DESIG_146(v) DESIG_145(v),[145]=v
#define DESIG__146(n,v) DESIG__145(n,v),[145]={DESIG_##n(v)}
#define EmptyPyGetSetDef_146 EmptyPyGetSetDef_145,[145]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_146 EmptyPyMethodDef_145,[145]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_146 EmptyPyMemberDef_145,[145]={NULL,0,0,0,NULL}
#define DESIG_147(v) DESIG_146(v),[146]=v
#define DESIG__147(n,v) DESIG__146(n,v),[146]={DESIG_##n(v)}
#define EmptyPyGetSetDef_147 EmptyPyGetSetDef_146,[146]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_147 EmptyPyMethodDef_146,[146]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_147 EmptyPyMemberDef_146,[146]={NULL,0,0,0,NULL}
#define DESIG_148(v) DESIG_147(v),[147]=v
#define DESIG__148(n,v) DESIG__147(n,v),[147]={DESIG_##n(v)}
#define EmptyPyGetSetDef_148 EmptyPyGetSetDef_147,[147]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_148 EmptyPyMethodDef_147,[147]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_148 EmptyPyMemberDef_147,[147]={NULL,0,0,0,NULL}
#define DESIG_149(v) DESIG_148(v),[148]=v
#define DESIG__149(n,v) DESIG__148(n,v),[148]={DESIG_##n(v)}
#define EmptyPyGetSetDef_149 EmptyPyGetSetDef_148,[148]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_149 EmptyPyMethodDef_148,[148]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_149 EmptyPyMemberDef_148,[148]={NULL,0,0,0,NULL}
#define DESIG_150(v) DESIG_149(v),[149]=v
#define DESIG__150(n,v) DESIG__149(n,v),[149]={DESIG_##n(v)}
#define EmptyPyGetSetDef_150 EmptyPyGetSetDef_149,[149]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_150 EmptyPyMethodDef_149,[149]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_150 EmptyPyMemberDef_149,[149]={NULL,0,0,0,NULL}
#define DESIG_151(v) DESIG_150(v),[150]=v
#define DESIG__151(n,v) DESIG__150(n,v),[150]={DESIG_##n(v)}
#define EmptyPyGetSetDef_151 EmptyPyGetSetDef_150,[150]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_151 EmptyPyMethodDef_150,[150]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_151 EmptyPyMemberDef_150,[150]={NULL,0,0,0,NULL}
#define DESIG_152(v) DESIG_151(v),[151]=v
#define DESIG__152(n,v) DESIG__151(n,v),[151]={DESIG_##n(v)}
#define EmptyPyGetSetDef_152 EmptyPyGetSetDef_151,[151]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_152 EmptyPyMethodDef_151,[151]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_152 EmptyPyMemberDef_151,[151]={NULL,0,0,0,NULL}
#define DESIG_153(v) DESIG_152(v),[152]=v
#define DESIG__153(n,v) DESIG__152(n,v),[152]={DESIG_##n(v)}
#define EmptyPyGetSetDef_153 EmptyPyGetSetDef_152,[152]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_153 EmptyPyMethodDef_152,[152]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_153 EmptyPyMemberDef_152,[152]={NULL,0,0,0,NULL}
#define DESIG_154(v) DESIG_153(v),[153]=v
#define DESIG__154(n,v) DESIG__153(n,v),[153]={DESIG_##n(v)}
#define EmptyPyGetSetDef_154 EmptyPyGetSetDef_153,[153]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_154 EmptyPyMethodDef_153,[153]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_154 EmptyPyMemberDef_153,[153]={NULL,0,0,0,NULL}
#define DESIG_155(v) DESIG_154(v),[154]=v
#define DESIG__155(n,v) DESIG__154(n,v),[154]={DESIG_##n(v)}
#define EmptyPyGetSetDef_155 EmptyPyGetSetDef_154,[154]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_155 EmptyPyMethodDef_154,[154]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_155 EmptyPyMemberDef_154,[154]={NULL,0,0,0,NULL}
#define DESIG_156(v) DESIG_155(v),[155]=v
#define DESIG__156(n,v) DESIG__155(n,v),[155]={DESIG_##n(v)}
#define EmptyPyGetSetDef_156 EmptyPyGetSetDef_155,[155]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_156 EmptyPyMethodDef_155,[155]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_156 EmptyPyMemberDef_155,[155]={NULL,0,0,0,NULL}
#define DESIG_157(v) DESIG_156(v),[156]=v
#define DESIG__157(n,v) DESIG__156(n,v),[156]={DESIG_##n(v)}
#define EmptyPyGetSetDef_157 EmptyPyGetSetDef_156,[156]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_157 EmptyPyMethodDef_156,[156]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_157 EmptyPyMemberDef_156,[156]={NULL,0,0,0,NULL}
#define DESIG_158(v) DESIG_157(v),[157]=v
#define DESIG__158(n,v) DESIG__157(n,v),[157]={DESIG_##n(v)}
#define EmptyPyGetSetDef_158 EmptyPyGetSetDef_157,[157]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_158 EmptyPyMethodDef_157,[157]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_158 EmptyPyMemberDef_157,[157]={NULL,0,0,0,NULL}
#define DESIG_159(v) DESIG_158(v),[158]=v
#define DESIG__159(n,v) DESIG__158(n,v),[158]={DESIG_##n(v)}
#define EmptyPyGetSetDef_159 EmptyPyGetSetDef_158,[158]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_159 EmptyPyMethodDef_158,[158]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_159 EmptyPyMemberDef_158,[158]={NULL,0,0,0,NULL}
#define DESIG_160(v) DESIG_159(v),[159]=v
#define DESIG__160(n,v) DESIG__159(n,v),[159]={DESIG_##n(v)}
#define EmptyPyGetSetDef_160 EmptyPyGetSetDef_159,[159]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_160 EmptyPyMethodDef_159,[159]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_160 EmptyPyMemberDef_159,[159]={NULL,0,0,0,NULL}
#define DESIG_161(v) DESIG_160(v),[160]=v
#define DESIG__161(n,v) DESIG__160(n,v),[160]={DESIG_##n(v)}
#define EmptyPyGetSetDef_161 EmptyPyGetSetDef_160,[160]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_161 EmptyPyMethodDef_160,[160]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_161 EmptyPyMemberDef_160,[160]={NULL,0,0,0,NULL}
#define DESIG_162(v) DESIG_161(v),[161]=v
#define DESIG__162(n,v) DESIG__161(n,v),[161]={DESIG_##n(v)}
#define EmptyPyGetSetDef_162 EmptyPyGetSetDef_161,[161]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_162 EmptyPyMethodDef_161,[161]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_162 EmptyPyMemberDef_161,[161]={NULL,0,0,0,NULL}
#define DESIG_163(v) DESIG_162(v),[162]=v
#define DESIG__163(n,v) DESIG__162(n,v),[162]={DESIG_##n(v)}
#define EmptyPyGetSetDef_163 EmptyPyGetSetDef_162,[162]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_163 EmptyPyMethodDef_162,[162]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_163 EmptyPyMemberDef_162,[162]={NULL,0,0,0,NULL}
#define DESIG_164(v) DESIG_163(v),[163]=v
#define DESIG__164(n,v) DESIG__163(n,v),[163]={DESIG_##n(v)}
#define EmptyPyGetSetDef_164 EmptyPyGetSetDef_163,[163]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_164 EmptyPyMethodDef_163,[163]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_164 EmptyPyMemberDef_163,[163]={NULL,0,0,0,NULL}
#define DESIG_165(v) DESIG_164(v),[164]=v
#define DESIG__165(n,v) DESIG__164(n,v),[164]={DESIG_##n(v)}
#define EmptyPyGetSetDef_165 EmptyPyGetSetDef_164,[164]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_165 EmptyPyMethodDef_164,[164]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_165 EmptyPyMemberDef_164,[164]={NULL,0,0,0,NULL}
#define DESIG_166(v) DESIG_165(v),[165]=v
#define DESIG__166(n,v) DESIG__165(n,v),[165]={DESIG_##n(v)}
#define EmptyPyGetSetDef_166 EmptyPyGetSetDef_165,[165]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_166 EmptyPyMethodDef_165,[165]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_166 EmptyPyMemberDef_165,[165]={NULL,0,0,0,NULL}
#define DESIG_167(v) DESIG_166(v),[166]=v
#define DESIG__167(n,v) DESIG__166(n,v),[166]={DESIG_##n(v)}
#define EmptyPyGetSetDef_167 EmptyPyGetSetDef_166,[166]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_167 EmptyPyMethodDef_166,[166]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_167 EmptyPyMemberDef_166,[166]={NULL,0,0,0,NULL}
#define DESIG_168(v) DESIG_167(v),[167]=v
#define DESIG__168(n,v) DESIG__167(n,v),[167]={DESIG_##n(v)}
#define EmptyPyGetSetDef_168 EmptyPyGetSetDef_167,[167]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_168 EmptyPyMethodDef_167,[167]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_168 EmptyPyMemberDef_167,[167]={NULL,0,0,0,NULL}
#define DESIG_169(v) DESIG_168(v),[168]=v
#define DESIG__169(n,v) DESIG__168(n,v),[168]={DESIG_##n(v)}
#define EmptyPyGetSetDef_169 EmptyPyGetSetDef_168,[168]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_169 EmptyPyMethodDef_168,[168]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_169 EmptyPyMemberDef_168,[168]={NULL,0,0,0,NULL}
#define DESIG_170(v) DESIG_169(v),[169]=v
#define DESIG__170(n,v) DESIG__169(n,v),[169]={DESIG_##n(v)}
#define EmptyPyGetSetDef_170 EmptyPyGetSetDef_169,[169]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_170 EmptyPyMethodDef_169,[169]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_170 EmptyPyMemberDef_169,[169]={NULL,0,0,0,NULL}
#define DESIG_171(v) DESIG_170(v),[170]=v
#define DESIG__171(n,v) DESIG__170(n,v),[170]={DESIG_##n(v)}
#define EmptyPyGetSetDef_171 EmptyPyGetSetDef_170,[170]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_171 EmptyPyMethodDef_170,[170]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_171 EmptyPyMemberDef_170,[170]={NULL,0,0,0,NULL}
#define DESIG_172(v) DESIG_171(v),[171]=v
#define DESIG__172(n,v) DESIG__171(n,v),[171]={DESIG_##n(v)}
#define EmptyPyGetSetDef_172 EmptyPyGetSetDef_171,[171]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_172 EmptyPyMethodDef_171,[171]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_172 EmptyPyMemberDef_171,[171]={NULL,0,0,0,NULL}
#define DESIG_173(v) DESIG_172(v),[172]=v
#define DESIG__173(n,v) DESIG__172(n,v),[172]={DESIG_##n(v)}
#define EmptyPyGetSetDef_173 EmptyPyGetSetDef_172,[172]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_173 EmptyPyMethodDef_172,[172]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_173 EmptyPyMemberDef_172,[172]={NULL,0,0,0,NULL}
#define DESIG_174(v) DESIG_173(v),[173]=v
#define DESIG__174(n,v) DESIG__173(n,v),[173]={DESIG_##n(v)}
#define EmptyPyGetSetDef_174 EmptyPyGetSetDef_173,[173]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_174 EmptyPyMethodDef_173,[173]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_174 EmptyPyMemberDef_173,[173]={NULL,0,0,0,NULL}
#define DESIG_175(v) DESIG_174(v),[174]=v
#define DESIG__175(n,v) DESIG__174(n,v),[174]={DESIG_##n(v)}
#define EmptyPyGetSetDef_175 EmptyPyGetSetDef_174,[174]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_175 EmptyPyMethodDef_174,[174]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_175 EmptyPyMemberDef_174,[174]={NULL,0,0,0,NULL}
#define DESIG_176(v) DESIG_175(v),[175]=v
#define DESIG__176(n,v) DESIG__175(n,v),[175]={DESIG_##n(v)}
#define EmptyPyGetSetDef_176 EmptyPyGetSetDef_175,[175]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_176 EmptyPyMethodDef_175,[175]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_176 EmptyPyMemberDef_175,[175]={NULL,0,0,0,NULL}
#define DESIG_177(v) DESIG_176(v),[176]=v
#define DESIG__177(n,v) DESIG__176(n,v),[176]={DESIG_##n(v)}
#define EmptyPyGetSetDef_177 EmptyPyGetSetDef_176,[176]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_177 EmptyPyMethodDef_176,[176]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_177 EmptyPyMemberDef_176,[176]={NULL,0,0,0,NULL}
#define DESIG_178(v) DESIG_177(v),[177]=v
#define DESIG__178(n,v) DESIG__177(n,v),[177]={DESIG_##n(v)}
#define EmptyPyGetSetDef_178 EmptyPyGetSetDef_177,[177]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_178 EmptyPyMethodDef_177,[177]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_178 EmptyPyMemberDef_177,[177]={NULL,0,0,0,NULL}
#define DESIG_179(v) DESIG_178(v),[178]=v
#define DESIG__179(n,v) DESIG__178(n,v),[178]={DESIG_##n(v)}
#define EmptyPyGetSetDef_179 EmptyPyGetSetDef_178,[178]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_179 EmptyPyMethodDef_178,[178]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_179 EmptyPyMemberDef_178,[178]={NULL,0,0,0,NULL}
#define DESIG_180(v) DESIG_179(v),[179]=v
#define DESIG__180(n,v) DESIG__179(n,v),[179]={DESIG_##n(v)}
#define EmptyPyGetSetDef_180 EmptyPyGetSetDef_179,[179]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_180 EmptyPyMethodDef_179,[179]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_180 EmptyPyMemberDef_179,[179]={NULL,0,0,0,NULL}
#define DESIG_181(v) DESIG_180(v),[180]=v
#define DESIG__181(n,v) DESIG__180(n,v),[180]={DESIG_##n(v)}
#define EmptyPyGetSetDef_181 EmptyPyGetSetDef_180,[180]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_181 EmptyPyMethodDef_180,[180]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_181 EmptyPyMemberDef_180,[180]={NULL,0,0,0,NULL}
#define DESIG_182(v) DESIG_181(v),[181]=v
#define DESIG__182(n,v) DESIG__181(n,v),[181]={DESIG_##n(v)}
#define EmptyPyGetSetDef_182 EmptyPyGetSetDef_181,[181]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_182 EmptyPyMethodDef_181,[181]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_182 EmptyPyMemberDef_181,[181]={NULL,0,0,0,NULL}
#define DESIG_183(v) DESIG_182(v),[182]=v
#define DESIG__183(n,v) DESIG__182(n,v),[182]={DESIG_##n(v)}
#define EmptyPyGetSetDef_183 EmptyPyGetSetDef_182,[182]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_183 EmptyPyMethodDef_182,[182]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_183 EmptyPyMemberDef_182,[182]={NULL,0,0,0,NULL}
#define DESIG_184(v) DESIG_183(v),[183]=v
#define DESIG__184(n,v) DESIG__183(n,v),[183]={DESIG_##n(v)}
#define EmptyPyGetSetDef_184 EmptyPyGetSetDef_183,[183]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_184 EmptyPyMethodDef_183,[183]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_184 EmptyPyMemberDef_183,[183]={NULL,0,0,0,NULL}
#define DESIG_185(v) DESIG_184(v),[184]=v
#define DESIG__185(n,v) DESIG__184(n,v),[184]={DESIG_##n(v)}
#define EmptyPyGetSetDef_185 EmptyPyGetSetDef_184,[184]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_185 EmptyPyMethodDef_184,[184]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_185 EmptyPyMemberDef_184,[184]={NULL,0,0,0,NULL}
#define DESIG_186(v) DESIG_185(v),[185]=v
#define DESIG__186(n,v) DESIG__185(n,v),[185]={DESIG_##n(v)}
#define EmptyPyGetSetDef_186 EmptyPyGetSetDef_185,[185]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_186 EmptyPyMethodDef_185,[185]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_186 EmptyPyMemberDef_185,[185]={NULL,0,0,0,NULL}
#define DESIG_187(v) DESIG_186(v),[186]=v
#define DESIG__187(n,v) DESIG__186(n,v),[186]={DESIG_##n(v)}
#define EmptyPyGetSetDef_187 EmptyPyGetSetDef_186,[186]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_187 EmptyPyMethodDef_186,[186]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_187 EmptyPyMemberDef_186,[186]={NULL,0,0,0,NULL}
#define DESIG_188(v) DESIG_187(v),[187]=v
#define DESIG__188(n,v) DESIG__187(n,v),[187]={DESIG_##n(v)}
#define EmptyPyGetSetDef_188 EmptyPyGetSetDef_187,[187]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_188 EmptyPyMethodDef_187,[187]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_188 EmptyPyMemberDef_187,[187]={NULL,0,0,0,NULL}
#define DESIG_189(v) DESIG_188(v),[188]=v
#define DESIG__189(n,v) DESIG__188(n,v),[188]={DESIG_##n(v)}
#define EmptyPyGetSetDef_189 EmptyPyGetSetDef_188,[188]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_189 EmptyPyMethodDef_188,[188]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_189 EmptyPyMemberDef_188,[188]={NULL,0,0,0,NULL}
#define DESIG_190(v) DESIG_189(v),[189]=v
#define DESIG__190(n,v) DESIG__189(n,v),[189]={DESIG_##n(v)}
#define EmptyPyGetSetDef_190 EmptyPyGetSetDef_189,[189]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_190 EmptyPyMethodDef_189,[189]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_190 EmptyPyMemberDef_189,[189]={NULL,0,0,0,NULL}
#define DESIG_191(v) DESIG_190(v),[190]=v
#define DESIG__191(n,v) DESIG__190(n,v),[190]={DESIG_##n(v)}
#define EmptyPyGetSetDef_191 EmptyPyGetSetDef_190,[190]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_191 EmptyPyMethodDef_190,[190]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_191 EmptyPyMemberDef_190,[190]={NULL,0,0,0,NULL}
#define DESIG_192(v) DESIG_191(v),[191]=v
#define DESIG__192(n,v) DESIG__191(n,v),[191]={DESIG_##n(v)}
#define EmptyPyGetSetDef_192 EmptyPyGetSetDef_191,[191]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_192 EmptyPyMethodDef_191,[191]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_192 EmptyPyMemberDef_191,[191]={NULL,0,0,0,NULL}
#define DESIG_193(v) DESIG_192(v),[192]=v
#define DESIG__193(n,v) DESIG__192(n,v),[192]={DESIG_##n(v)}
#define EmptyPyGetSetDef_193 EmptyPyGetSetDef_192,[192]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_193 EmptyPyMethodDef_192,[192]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_193 EmptyPyMemberDef_192,[192]={NULL,0,0,0,NULL}
#define DESIG_194(v) DESIG_193(v),[193]=v
#define DESIG__194(n,v) DESIG__193(n,v),[193]={DESIG_##n(v)}
#define EmptyPyGetSetDef_194 EmptyPyGetSetDef_193,[193]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_194 EmptyPyMethodDef_193,[193]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_194 EmptyPyMemberDef_193,[193]={NULL,0,0,0,NULL}
#define DESIG_195(v) DESIG_194(v),[194]=v
#define DESIG__195(n,v) DESIG__194(n,v),[194]={DESIG_##n(v)}
#define EmptyPyGetSetDef_195 EmptyPyGetSetDef_194,[194]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_195 EmptyPyMethodDef_194,[194]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_195 EmptyPyMemberDef_194,[194]={NULL,0,0,0,NULL}
#define DESIG_196(v) DESIG_195(v),[195]=v
#define DESIG__196(n,v) DESIG__195(n,v),[195]={DESIG_##n(v)}
#define EmptyPyGetSetDef_196 EmptyPyGetSetDef_195,[195]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_196 EmptyPyMethodDef_195,[195]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_196 EmptyPyMemberDef_195,[195]={NULL,0,0,0,NULL}
#define DESIG_197(v) DESIG_196(v),[196]=v
#define DESIG__197(n,v) DESIG__196(n,v),[196]={DESIG_##n(v)}
#define EmptyPyGetSetDef_197 EmptyPyGetSetDef_196,[196]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_197 EmptyPyMethodDef_196,[196]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_197 EmptyPyMemberDef_196,[196]={NULL,0,0,0,NULL}
#define DESIG_198(v) DESIG_197(v),[197]=v
#define DESIG__198(n,v) DESIG__197(n,v),[197]={DESIG_##n(v)}
#define EmptyPyGetSetDef_198 EmptyPyGetSetDef_197,[197]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_198 EmptyPyMethodDef_197,[197]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_198 EmptyPyMemberDef_197,[197]={NULL,0,0,0,NULL}
#define DESIG_199(v) DESIG_198(v),[198]=v
#define DESIG__199(n,v) DESIG__198(n,v),[198]={DESIG_##n(v)}
#define EmptyPyGetSetDef_199 EmptyPyGetSetDef_198,[198]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_199 EmptyPyMethodDef_198,[198]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_199 EmptyPyMemberDef_198,[198]={NULL,0,0,0,NULL}
#define DESIG_200(v) DESIG_199(v),[199]=v
#define DESIG__200(n,v) DESIG__199(n,v),[199]={DESIG_##n(v)}
#define EmptyPyGetSetDef_200 EmptyPyGetSetDef_199,[199]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_200 EmptyPyMethodDef_199,[199]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_200 EmptyPyMemberDef_199,[199]={NULL,0,0,0,NULL}
#define DESIG_201(v) DESIG_200(v),[200]=v
#define DESIG__201(n,v) DESIG__200(n,v),[200]={DESIG_##n(v)}
#define EmptyPyGetSetDef_201 EmptyPyGetSetDef_200,[200]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_201 EmptyPyMethodDef_200,[200]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_201 EmptyPyMemberDef_200,[200]={NULL,0,0,0,NULL}
#define DESIG_202(v) DESIG_201(v),[201]=v
#define DESIG__202(n,v) DESIG__201(n,v),[201]={DESIG_##n(v)}
#define EmptyPyGetSetDef_202 EmptyPyGetSetDef_201,[201]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_202 EmptyPyMethodDef_201,[201]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_202 EmptyPyMemberDef_201,[201]={NULL,0,0,0,NULL}
#define DESIG_203(v) DESIG_202(v),[202]=v
#define DESIG__203(n,v) DESIG__202(n,v),[202]={DESIG_##n(v)}
#define EmptyPyGetSetDef_203 EmptyPyGetSetDef_202,[202]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_203 EmptyPyMethodDef_202,[202]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_203 EmptyPyMemberDef_202,[202]={NULL,0,0,0,NULL}
#define DESIG_204(v) DESIG_203(v),[203]=v
#define DESIG__204(n,v) DESIG__203(n,v),[203]={DESIG_##n(v)}
#define EmptyPyGetSetDef_204 EmptyPyGetSetDef_203,[203]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_204 EmptyPyMethodDef_203,[203]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_204 EmptyPyMemberDef_203,[203]={NULL,0,0,0,NULL}
#define DESIG_205(v) DESIG_204(v),[204]=v
#define DESIG__205(n,v) DESIG__204(n,v),[204]={DESIG_##n(v)}
#define EmptyPyGetSetDef_205 EmptyPyGetSetDef_204,[204]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_205 EmptyPyMethodDef_204,[204]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_205 EmptyPyMemberDef_204,[204]={NULL,0,0,0,NULL}
#define DESIG_206(v) DESIG_205(v),[205]=v
#define DESIG__206(n,v) DESIG__205(n,v),[205]={DESIG_##n(v)}
#define EmptyPyGetSetDef_206 EmptyPyGetSetDef_205,[205]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_206 EmptyPyMethodDef_205,[205]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_206 EmptyPyMemberDef_205,[205]={NULL,0,0,0,NULL}
#define DESIG_207(v) DESIG_206(v),[206]=v
#define DESIG__207(n,v) DESIG__206(n,v),[206]={DESIG_##n(v)}
#define EmptyPyGetSetDef_207 EmptyPyGetSetDef_206,[206]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_207 EmptyPyMethodDef_206,[206]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_207 EmptyPyMemberDef_206,[206]={NULL,0,0,0,NULL}
#define DESIG_208(v) DESIG_207(v),[207]=v
#define DESIG__208(n,v) DESIG__207(n,v),[207]={DESIG_##n(v)}
#define EmptyPyGetSetDef_208 EmptyPyGetSetDef_207,[207]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_208 EmptyPyMethodDef_207,[207]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_208 EmptyPyMemberDef_207,[207]={NULL,0,0,0,NULL}
#define DESIG_209(v) DESIG_208(v),[208]=v
#define DESIG__209(n,v) DESIG__208(n,v),[208]={DESIG_##n(v)}
#define EmptyPyGetSetDef_209 EmptyPyGetSetDef_208,[208]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_209 EmptyPyMethodDef_208,[208]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_209 EmptyPyMemberDef_208,[208]={NULL,0,0,0,NULL}
#define DESIG_210(v) DESIG_209(v),[209]=v
#define DESIG__210(n,v) DESIG__209(n,v),[209]={DESIG_##n(v)}
#define EmptyPyGetSetDef_210 EmptyPyGetSetDef_209,[209]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_210 EmptyPyMethodDef_209,[209]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_210 EmptyPyMemberDef_209,[209]={NULL,0,0,0,NULL}
#define DESIG_211(v) DESIG_210(v),[210]=v
#define DESIG__211(n,v) DESIG__210(n,v),[210]={DESIG_##n(v)}
#define EmptyPyGetSetDef_211 EmptyPyGetSetDef_210,[210]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_211 EmptyPyMethodDef_210,[210]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_211 EmptyPyMemberDef_210,[210]={NULL,0,0,0,NULL}
#define DESIG_212(v) DESIG_211(v),[211]=v
#define DESIG__212(n,v) DESIG__211(n,v),[211]={DESIG_##n(v)}
#define EmptyPyGetSetDef_212 EmptyPyGetSetDef_211,[211]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_212 EmptyPyMethodDef_211,[211]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_212 EmptyPyMemberDef_211,[211]={NULL,0,0,0,NULL}
#define DESIG_213(v) DESIG_212(v),[212]=v
#define DESIG__213(n,v) DESIG__212(n,v),[212]={DESIG_##n(v)}
#define EmptyPyGetSetDef_213 EmptyPyGetSetDef_212,[212]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_213 EmptyPyMethodDef_212,[212]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_213 EmptyPyMemberDef_212,[212]={NULL,0,0,0,NULL}
#define DESIG_214(v) DESIG_213(v),[213]=v
#define DESIG__214(n,v) DESIG__213(n,v),[213]={DESIG_##n(v)}
#define EmptyPyGetSetDef_214 EmptyPyGetSetDef_213,[213]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_214 EmptyPyMethodDef_213,[213]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_214 EmptyPyMemberDef_213,[213]={NULL,0,0,0,NULL}
#define DESIG_215(v) DESIG_214(v),[214]=v
#define DESIG__215(n,v) DESIG__214(n,v),[214]={DESIG_##n(v)}
#define EmptyPyGetSetDef_215 EmptyPyGetSetDef_214,[214]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_215 EmptyPyMethodDef_214,[214]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_215 EmptyPyMemberDef_214,[214]={NULL,0,0,0,NULL}
#define DESIG_216(v) DESIG_215(v),[215]=v
#define DESIG__216(n,v) DESIG__215(n,v),[215]={DESIG_##n(v)}
#define EmptyPyGetSetDef_216 EmptyPyGetSetDef_215,[215]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_216 EmptyPyMethodDef_215,[215]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_216 EmptyPyMemberDef_215,[215]={NULL,0,0,0,NULL}
#define DESIG_217(v) DESIG_216(v),[216]=v
#define DESIG__217(n,v) DESIG__216(n,v),[216]={DESIG_##n(v)}
#define EmptyPyGetSetDef_217 EmptyPyGetSetDef_216,[216]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_217 EmptyPyMethodDef_216,[216]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_217 EmptyPyMemberDef_216,[216]={NULL,0,0,0,NULL}
#define DESIG_218(v) DESIG_217(v),[217]=v
#define DESIG__218(n,v) DESIG__217(n,v),[217]={DESIG_##n(v)}
#define EmptyPyGetSetDef_218 EmptyPyGetSetDef_217,[217]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_218 EmptyPyMethodDef_217,[217]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_218 EmptyPyMemberDef_217,[217]={NULL,0,0,0,NULL}
#define DESIG_219(v) DESIG_218(v),[218]=v
#define DESIG__219(n,v) DESIG__218(n,v),[218]={DESIG_##n(v)}
#define EmptyPyGetSetDef_219 EmptyPyGetSetDef_218,[218]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_219 EmptyPyMethodDef_218,[218]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_219 EmptyPyMemberDef_218,[218]={NULL,0,0,0,NULL}
#define DESIG_220(v) DESIG_219(v),[219]=v
#define DESIG__220(n,v) DESIG__219(n,v),[219]={DESIG_##n(v)}
#define EmptyPyGetSetDef_220 EmptyPyGetSetDef_219,[219]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_220 EmptyPyMethodDef_219,[219]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_220 EmptyPyMemberDef_219,[219]={NULL,0,0,0,NULL}
#define DESIG_221(v) DESIG_220(v),[220]=v
#define DESIG__221(n,v) DESIG__220(n,v),[220]={DESIG_##n(v)}
#define EmptyPyGetSetDef_221 EmptyPyGetSetDef_220,[220]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_221 EmptyPyMethodDef_220,[220]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_221 EmptyPyMemberDef_220,[220]={NULL,0,0,0,NULL}
#define DESIG_222(v) DESIG_221(v),[221]=v
#define DESIG__222(n,v) DESIG__221(n,v),[221]={DESIG_##n(v)}
#define EmptyPyGetSetDef_222 EmptyPyGetSetDef_221,[221]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_222 EmptyPyMethodDef_221,[221]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_222 EmptyPyMemberDef_221,[221]={NULL,0,0,0,NULL}
#define DESIG_223(v) DESIG_222(v),[222]=v
#define DESIG__223(n,v) DESIG__222(n,v),[222]={DESIG_##n(v)}
#define EmptyPyGetSetDef_223 EmptyPyGetSetDef_222,[222]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_223 EmptyPyMethodDef_222,[222]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_223 EmptyPyMemberDef_222,[222]={NULL,0,0,0,NULL}
#define DESIG_224(v) DESIG_223(v),[223]=v
#define DESIG__224(n,v) DESIG__223(n,v),[223]={DESIG_##n(v)}
#define EmptyPyGetSetDef_224 EmptyPyGetSetDef_223,[223]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_224 EmptyPyMethodDef_223,[223]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_224 EmptyPyMemberDef_223,[223]={NULL,0,0,0,NULL}
#define DESIG_225(v) DESIG_224(v),[224]=v
#define DESIG__225(n,v) DESIG__224(n,v),[224]={DESIG_##n(v)}
#define EmptyPyGetSetDef_225 EmptyPyGetSetDef_224,[224]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_225 EmptyPyMethodDef_224,[224]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_225 EmptyPyMemberDef_224,[224]={NULL,0,0,0,NULL}
#define DESIG_226(v) DESIG_225(v),[225]=v
#define DESIG__226(n,v) DESIG__225(n,v),[225]={DESIG_##n(v)}
#define EmptyPyGetSetDef_226 EmptyPyGetSetDef_225,[225]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_226 EmptyPyMethodDef_225,[225]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_226 EmptyPyMemberDef_225,[225]={NULL,0,0,0,NULL}
#define DESIG_227(v) DESIG_226(v),[226]=v
#define DESIG__227(n,v) DESIG__226(n,v),[226]={DESIG_##n(v)}
#define EmptyPyGetSetDef_227 EmptyPyGetSetDef_226,[226]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_227 EmptyPyMethodDef_226,[226]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_227 EmptyPyMemberDef_226,[226]={NULL,0,0,0,NULL}
#define DESIG_228(v) DESIG_227(v),[227]=v
#define DESIG__228(n,v) DESIG__227(n,v),[227]={DESIG_##n(v)}
#define EmptyPyGetSetDef_228 EmptyPyGetSetDef_227,[227]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_228 EmptyPyMethodDef_227,[227]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_228 EmptyPyMemberDef_227,[227]={NULL,0,0,0,NULL}
#define DESIG_229(v) DESIG_228(v),[228]=v
#define DESIG__229(n,v) DESIG__228(n,v),[228]={DESIG_##n(v)}
#define EmptyPyGetSetDef_229 EmptyPyGetSetDef_228,[228]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_229 EmptyPyMethodDef_228,[228]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_229 EmptyPyMemberDef_228,[228]={NULL,0,0,0,NULL}
#define DESIG_230(v) DESIG_229(v),[229]=v
#define DESIG__230(n,v) DESIG__229(n,v),[229]={DESIG_##n(v)}
#define EmptyPyGetSetDef_230 EmptyPyGetSetDef_229,[229]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_230 EmptyPyMethodDef_229,[229]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_230 EmptyPyMemberDef_229,[229]={NULL,0,0,0,NULL}
#define DESIG_231(v) DESIG_230(v),[230]=v
#define DESIG__231(n,v) DESIG__230(n,v),[230]={DESIG_##n(v)}
#define EmptyPyGetSetDef_231 EmptyPyGetSetDef_230,[230]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_231 EmptyPyMethodDef_230,[230]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_231 EmptyPyMemberDef_230,[230]={NULL,0,0,0,NULL}
#define DESIG_232(v) DESIG_231(v),[231]=v
#define DESIG__232(n,v) DESIG__231(n,v),[231]={DESIG_##n(v)}
#define EmptyPyGetSetDef_232 EmptyPyGetSetDef_231,[231]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_232 EmptyPyMethodDef_231,[231]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_232 EmptyPyMemberDef_231,[231]={NULL,0,0,0,NULL}
#define DESIG_233(v) DESIG_232(v),[232]=v
#define DESIG__233(n,v) DESIG__232(n,v),[232]={DESIG_##n(v)}
#define EmptyPyGetSetDef_233 EmptyPyGetSetDef_232,[232]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_233 EmptyPyMethodDef_232,[232]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_233 EmptyPyMemberDef_232,[232]={NULL,0,0,0,NULL}
#define DESIG_234(v) DESIG_233(v),[233]=v
#define DESIG__234(n,v) DESIG__233(n,v),[233]={DESIG_##n(v)}
#define EmptyPyGetSetDef_234 EmptyPyGetSetDef_233,[233]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_234 EmptyPyMethodDef_233,[233]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_234 EmptyPyMemberDef_233,[233]={NULL,0,0,0,NULL}
#define DESIG_235(v) DESIG_234(v),[234]=v
#define DESIG__235(n,v) DESIG__234(n,v),[234]={DESIG_##n(v)}
#define EmptyPyGetSetDef_235 EmptyPyGetSetDef_234,[234]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_235 EmptyPyMethodDef_234,[234]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_235 EmptyPyMemberDef_234,[234]={NULL,0,0,0,NULL}
#define DESIG_236(v) DESIG_235(v),[235]=v
#define DESIG__236(n,v) DESIG__235(n,v),[235]={DESIG_##n(v)}
#define EmptyPyGetSetDef_236 EmptyPyGetSetDef_235,[235]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_236 EmptyPyMethodDef_235,[235]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_236 EmptyPyMemberDef_235,[235]={NULL,0,0,0,NULL}
#define DESIG_237(v) DESIG_236(v),[236]=v
#define DESIG__237(n,v) DESIG__236(n,v),[236]={DESIG_##n(v)}
#define EmptyPyGetSetDef_237 EmptyPyGetSetDef_236,[236]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_237 EmptyPyMethodDef_236,[236]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_237 EmptyPyMemberDef_236,[236]={NULL,0,0,0,NULL}
#define DESIG_238(v) DESIG_237(v),[237]=v
#define DESIG__238(n,v) DESIG__237(n,v),[237]={DESIG_##n(v)}
#define EmptyPyGetSetDef_238 EmptyPyGetSetDef_237,[237]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_238 EmptyPyMethodDef_237,[237]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_238 EmptyPyMemberDef_237,[237]={NULL,0,0,0,NULL}
#define DESIG_239(v) DESIG_238(v),[238]=v
#define DESIG__239(n,v) DESIG__238(n,v),[238]={DESIG_##n(v)}
#define EmptyPyGetSetDef_239 EmptyPyGetSetDef_238,[238]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_239 EmptyPyMethodDef_238,[238]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_239 EmptyPyMemberDef_238,[238]={NULL,0,0,0,NULL}
#define DESIG_240(v) DESIG_239(v),[239]=v
#define DESIG__240(n,v) DESIG__239(n,v),[239]={DESIG_##n(v)}
#define EmptyPyGetSetDef_240 EmptyPyGetSetDef_239,[239]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_240 EmptyPyMethodDef_239,[239]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_240 EmptyPyMemberDef_239,[239]={NULL,0,0,0,NULL}
#define DESIG_241(v) DESIG_240(v),[240]=v
#define DESIG__241(n,v) DESIG__240(n,v),[240]={DESIG_##n(v)}
#define EmptyPyGetSetDef_241 EmptyPyGetSetDef_240,[240]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_241 EmptyPyMethodDef_240,[240]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_241 EmptyPyMemberDef_240,[240]={NULL,0,0,0,NULL}
#define DESIG_242(v) DESIG_241(v),[241]=v
#define DESIG__242(n,v) DESIG__241(n,v),[241]={DESIG_##n(v)}
#define EmptyPyGetSetDef_242 EmptyPyGetSetDef_241,[241]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_242 EmptyPyMethodDef_241,[241]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_242 EmptyPyMemberDef_241,[241]={NULL,0,0,0,NULL}
#define DESIG_243(v) DESIG_242(v),[242]=v
#define DESIG__243(n,v) DESIG__242(n,v),[242]={DESIG_##n(v)}
#define EmptyPyGetSetDef_243 EmptyPyGetSetDef_242,[242]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_243 EmptyPyMethodDef_242,[242]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_243 EmptyPyMemberDef_242,[242]={NULL,0,0,0,NULL}
#define DESIG_244(v) DESIG_243(v),[243]=v
#define DESIG__244(n,v) DESIG__243(n,v),[243]={DESIG_##n(v)}
#define EmptyPyGetSetDef_244 EmptyPyGetSetDef_243,[243]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_244 EmptyPyMethodDef_243,[243]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_244 EmptyPyMemberDef_243,[243]={NULL,0,0,0,NULL}
#define DESIG_245(v) DESIG_244(v),[244]=v
#define DESIG__245(n,v) DESIG__244(n,v),[244]={DESIG_##n(v)}
#define EmptyPyGetSetDef_245 EmptyPyGetSetDef_244,[244]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_245 EmptyPyMethodDef_244,[244]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_245 EmptyPyMemberDef_244,[244]={NULL,0,0,0,NULL}
#define DESIG_246(v) DESIG_245(v),[245]=v
#define DESIG__246(n,v) DESIG__245(n,v),[245]={DESIG_##n(v)}
#define EmptyPyGetSetDef_246 EmptyPyGetSetDef_245,[245]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_246 EmptyPyMethodDef_245,[245]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_246 EmptyPyMemberDef_245,[245]={NULL,0,0,0,NULL}
#define DESIG_247(v) DESIG_246(v),[246]=v
#define DESIG__247(n,v) DESIG__246(n,v),[246]={DESIG_##n(v)}
#define EmptyPyGetSetDef_247 EmptyPyGetSetDef_246,[246]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_247 EmptyPyMethodDef_246,[246]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_247 EmptyPyMemberDef_246,[246]={NULL,0,0,0,NULL}
#define DESIG_248(v) DESIG_247(v),[247]=v
#define DESIG__248(n,v) DESIG__247(n,v),[247]={DESIG_##n(v)}
#define EmptyPyGetSetDef_248 EmptyPyGetSetDef_247,[247]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_248 EmptyPyMethodDef_247,[247]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_248 EmptyPyMemberDef_247,[247]={NULL,0,0,0,NULL}
#define DESIG_249(v) DESIG_248(v),[248]=v
#define DESIG__249(n,v) DESIG__248(n,v),[248]={DESIG_##n(v)}
#define EmptyPyGetSetDef_249 EmptyPyGetSetDef_248,[248]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_249 EmptyPyMethodDef_248,[248]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_249 EmptyPyMemberDef_248,[248]={NULL,0,0,0,NULL}
#define DESIG_250(v) DESIG_249(v),[249]=v
#define DESIG__250(n,v) DESIG__249(n,v),[249]={DESIG_##n(v)}
#define EmptyPyGetSetDef_250 EmptyPyGetSetDef_249,[249]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_250 EmptyPyMethodDef_249,[249]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_250 EmptyPyMemberDef_249,[249]={NULL,0,0,0,NULL}
#define DESIG_251(v) DESIG_250(v),[250]=v
#define DESIG__251(n,v) DESIG__250(n,v),[250]={DESIG_##n(v)}
#define EmptyPyGetSetDef_251 EmptyPyGetSetDef_250,[250]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_251 EmptyPyMethodDef_250,[250]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_251 EmptyPyMemberDef_250,[250]={NULL,0,0,0,NULL}
#define DESIG_252(v) DESIG_251(v),[251]=v
#define DESIG__252(n,v) DESIG__251(n,v),[251]={DESIG_##n(v)}
#define EmptyPyGetSetDef_252 EmptyPyGetSetDef_251,[251]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_252 EmptyPyMethodDef_251,[251]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_252 EmptyPyMemberDef_251,[251]={NULL,0,0,0,NULL}
#define DESIG_253(v) DESIG_252(v),[252]=v
#define DESIG__253(n,v) DESIG__252(n,v),[252]={DESIG_##n(v)}
#define EmptyPyGetSetDef_253 EmptyPyGetSetDef_252,[252]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_253 EmptyPyMethodDef_252,[252]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_253 EmptyPyMemberDef_252,[252]={NULL,0,0,0,NULL}
#define DESIG_254(v) DESIG_253(v),[253]=v
#define DESIG__254(n,v) DESIG__253(n,v),[253]={DESIG_##n(v)}
#define EmptyPyGetSetDef_254 EmptyPyGetSetDef_253,[253]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_254 EmptyPyMethodDef_253,[253]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_254 EmptyPyMemberDef_253,[253]={NULL,0,0,0,NULL}
#define DESIG_255(v) DESIG_254(v),[254]=v
#define DESIG__255(n,v) DESIG__254(n,v),[254]={DESIG_##n(v)}
#define EmptyPyGetSetDef_255 EmptyPyGetSetDef_254,[254]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_255 EmptyPyMethodDef_254,[254]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_255 EmptyPyMemberDef_254,[254]={NULL,0,0,0,NULL}
#define DESIG_256(v) DESIG_255(v),[255]=v
#define DESIG__256(n,v) DESIG__255(n,v),[255]={DESIG_##n(v)}
#define EmptyPyGetSetDef_256 EmptyPyGetSetDef_255,[255]={NULL,NULL,NULL,NULL,NULL}
#define EmptyPyMethodDef_256 EmptyPyMethodDef_255,[255]={NULL,NULL,0,NULL}
#define EmptyPyMemberDef_256 EmptyPyMemberDef_255,[255]={NULL,0,0,0,NULL}
#define DESIG(m,v) DESIG_##m(v)
#define DESIG2(m,n,v) DESIG__##m(n,v)
#define EmptyPyGetSetDef(m) {EmptyPyGetSetDef_##m}
#define EmptyPyMethodDef(m) {EmptyPyMethodDef_##m}
#define EmptyPyMemberDef(m) {EmptyPyMemberDef_##m}



#define EmptyModule {PyModuleDef_HEAD_INIT,NULL,NULL,-1,NULL,NULL,NULL,NULL}


/*
//the code that was used to generate the following bullshit macros:
 
#include<stdio.h>
int main()
{
 
    int max_n=256;
    for(int i=0;i<max_n+1;i++)
    {
        if(i==0) printf("#define POTFIT_THRMO_0 2,\"ERR\",err\n");
        else printf("#define POTFIT_THRMO_%d POTFIT_THRMO_%d,names[%d],errs[%d]\n",i,i-1,i-1,i-1);
 
    }
 
    for(int i=0;i<max_n+1;i++)
    {
        if(i!=max_n) printf("#define RET_POTFIT_THRMO_%d if(nconfigs==%d) return ThermoDynamics(POTFIT_THRMO_%d); RET_POTFIT_THRMO_%d\n",i,i,i,i+1);
        else printf("#define RET_POTFIT_THRMO_%d return ThermoDynamics(POTFIT_THRMO_%d);\n",i,i);
    }
 
    printf("#define RET_POTFIT_THRMO RET_POTFIT_THRMO_0\n");
 
    return 0;
}
 */









#define POTFIT_THRMO_0 2,"ERR",err
#define POTFIT_THRMO_1 POTFIT_THRMO_0,names[0],errs[0]
#define POTFIT_THRMO_2 POTFIT_THRMO_1,names[1],errs[1]
#define POTFIT_THRMO_3 POTFIT_THRMO_2,names[2],errs[2]
#define POTFIT_THRMO_4 POTFIT_THRMO_3,names[3],errs[3]
#define POTFIT_THRMO_5 POTFIT_THRMO_4,names[4],errs[4]
#define POTFIT_THRMO_6 POTFIT_THRMO_5,names[5],errs[5]
#define POTFIT_THRMO_7 POTFIT_THRMO_6,names[6],errs[6]
#define POTFIT_THRMO_8 POTFIT_THRMO_7,names[7],errs[7]
#define POTFIT_THRMO_9 POTFIT_THRMO_8,names[8],errs[8]
#define POTFIT_THRMO_10 POTFIT_THRMO_9,names[9],errs[9]
#define POTFIT_THRMO_11 POTFIT_THRMO_10,names[10],errs[10]
#define POTFIT_THRMO_12 POTFIT_THRMO_11,names[11],errs[11]
#define POTFIT_THRMO_13 POTFIT_THRMO_12,names[12],errs[12]
#define POTFIT_THRMO_14 POTFIT_THRMO_13,names[13],errs[13]
#define POTFIT_THRMO_15 POTFIT_THRMO_14,names[14],errs[14]
#define POTFIT_THRMO_16 POTFIT_THRMO_15,names[15],errs[15]
#define POTFIT_THRMO_17 POTFIT_THRMO_16,names[16],errs[16]
#define POTFIT_THRMO_18 POTFIT_THRMO_17,names[17],errs[17]
#define POTFIT_THRMO_19 POTFIT_THRMO_18,names[18],errs[18]
#define POTFIT_THRMO_20 POTFIT_THRMO_19,names[19],errs[19]
#define POTFIT_THRMO_21 POTFIT_THRMO_20,names[20],errs[20]
#define POTFIT_THRMO_22 POTFIT_THRMO_21,names[21],errs[21]
#define POTFIT_THRMO_23 POTFIT_THRMO_22,names[22],errs[22]
#define POTFIT_THRMO_24 POTFIT_THRMO_23,names[23],errs[23]
#define POTFIT_THRMO_25 POTFIT_THRMO_24,names[24],errs[24]
#define POTFIT_THRMO_26 POTFIT_THRMO_25,names[25],errs[25]
#define POTFIT_THRMO_27 POTFIT_THRMO_26,names[26],errs[26]
#define POTFIT_THRMO_28 POTFIT_THRMO_27,names[27],errs[27]
#define POTFIT_THRMO_29 POTFIT_THRMO_28,names[28],errs[28]
#define POTFIT_THRMO_30 POTFIT_THRMO_29,names[29],errs[29]
#define POTFIT_THRMO_31 POTFIT_THRMO_30,names[30],errs[30]
#define POTFIT_THRMO_32 POTFIT_THRMO_31,names[31],errs[31]
#define POTFIT_THRMO_33 POTFIT_THRMO_32,names[32],errs[32]
#define POTFIT_THRMO_34 POTFIT_THRMO_33,names[33],errs[33]
#define POTFIT_THRMO_35 POTFIT_THRMO_34,names[34],errs[34]
#define POTFIT_THRMO_36 POTFIT_THRMO_35,names[35],errs[35]
#define POTFIT_THRMO_37 POTFIT_THRMO_36,names[36],errs[36]
#define POTFIT_THRMO_38 POTFIT_THRMO_37,names[37],errs[37]
#define POTFIT_THRMO_39 POTFIT_THRMO_38,names[38],errs[38]
#define POTFIT_THRMO_40 POTFIT_THRMO_39,names[39],errs[39]
#define POTFIT_THRMO_41 POTFIT_THRMO_40,names[40],errs[40]
#define POTFIT_THRMO_42 POTFIT_THRMO_41,names[41],errs[41]
#define POTFIT_THRMO_43 POTFIT_THRMO_42,names[42],errs[42]
#define POTFIT_THRMO_44 POTFIT_THRMO_43,names[43],errs[43]
#define POTFIT_THRMO_45 POTFIT_THRMO_44,names[44],errs[44]
#define POTFIT_THRMO_46 POTFIT_THRMO_45,names[45],errs[45]
#define POTFIT_THRMO_47 POTFIT_THRMO_46,names[46],errs[46]
#define POTFIT_THRMO_48 POTFIT_THRMO_47,names[47],errs[47]
#define POTFIT_THRMO_49 POTFIT_THRMO_48,names[48],errs[48]
#define POTFIT_THRMO_50 POTFIT_THRMO_49,names[49],errs[49]
#define POTFIT_THRMO_51 POTFIT_THRMO_50,names[50],errs[50]
#define POTFIT_THRMO_52 POTFIT_THRMO_51,names[51],errs[51]
#define POTFIT_THRMO_53 POTFIT_THRMO_52,names[52],errs[52]
#define POTFIT_THRMO_54 POTFIT_THRMO_53,names[53],errs[53]
#define POTFIT_THRMO_55 POTFIT_THRMO_54,names[54],errs[54]
#define POTFIT_THRMO_56 POTFIT_THRMO_55,names[55],errs[55]
#define POTFIT_THRMO_57 POTFIT_THRMO_56,names[56],errs[56]
#define POTFIT_THRMO_58 POTFIT_THRMO_57,names[57],errs[57]
#define POTFIT_THRMO_59 POTFIT_THRMO_58,names[58],errs[58]
#define POTFIT_THRMO_60 POTFIT_THRMO_59,names[59],errs[59]
#define POTFIT_THRMO_61 POTFIT_THRMO_60,names[60],errs[60]
#define POTFIT_THRMO_62 POTFIT_THRMO_61,names[61],errs[61]
#define POTFIT_THRMO_63 POTFIT_THRMO_62,names[62],errs[62]
#define POTFIT_THRMO_64 POTFIT_THRMO_63,names[63],errs[63]
#define POTFIT_THRMO_65 POTFIT_THRMO_64,names[64],errs[64]
#define POTFIT_THRMO_66 POTFIT_THRMO_65,names[65],errs[65]
#define POTFIT_THRMO_67 POTFIT_THRMO_66,names[66],errs[66]
#define POTFIT_THRMO_68 POTFIT_THRMO_67,names[67],errs[67]
#define POTFIT_THRMO_69 POTFIT_THRMO_68,names[68],errs[68]
#define POTFIT_THRMO_70 POTFIT_THRMO_69,names[69],errs[69]
#define POTFIT_THRMO_71 POTFIT_THRMO_70,names[70],errs[70]
#define POTFIT_THRMO_72 POTFIT_THRMO_71,names[71],errs[71]
#define POTFIT_THRMO_73 POTFIT_THRMO_72,names[72],errs[72]
#define POTFIT_THRMO_74 POTFIT_THRMO_73,names[73],errs[73]
#define POTFIT_THRMO_75 POTFIT_THRMO_74,names[74],errs[74]
#define POTFIT_THRMO_76 POTFIT_THRMO_75,names[75],errs[75]
#define POTFIT_THRMO_77 POTFIT_THRMO_76,names[76],errs[76]
#define POTFIT_THRMO_78 POTFIT_THRMO_77,names[77],errs[77]
#define POTFIT_THRMO_79 POTFIT_THRMO_78,names[78],errs[78]
#define POTFIT_THRMO_80 POTFIT_THRMO_79,names[79],errs[79]
#define POTFIT_THRMO_81 POTFIT_THRMO_80,names[80],errs[80]
#define POTFIT_THRMO_82 POTFIT_THRMO_81,names[81],errs[81]
#define POTFIT_THRMO_83 POTFIT_THRMO_82,names[82],errs[82]
#define POTFIT_THRMO_84 POTFIT_THRMO_83,names[83],errs[83]
#define POTFIT_THRMO_85 POTFIT_THRMO_84,names[84],errs[84]
#define POTFIT_THRMO_86 POTFIT_THRMO_85,names[85],errs[85]
#define POTFIT_THRMO_87 POTFIT_THRMO_86,names[86],errs[86]
#define POTFIT_THRMO_88 POTFIT_THRMO_87,names[87],errs[87]
#define POTFIT_THRMO_89 POTFIT_THRMO_88,names[88],errs[88]
#define POTFIT_THRMO_90 POTFIT_THRMO_89,names[89],errs[89]
#define POTFIT_THRMO_91 POTFIT_THRMO_90,names[90],errs[90]
#define POTFIT_THRMO_92 POTFIT_THRMO_91,names[91],errs[91]
#define POTFIT_THRMO_93 POTFIT_THRMO_92,names[92],errs[92]
#define POTFIT_THRMO_94 POTFIT_THRMO_93,names[93],errs[93]
#define POTFIT_THRMO_95 POTFIT_THRMO_94,names[94],errs[94]
#define POTFIT_THRMO_96 POTFIT_THRMO_95,names[95],errs[95]
#define POTFIT_THRMO_97 POTFIT_THRMO_96,names[96],errs[96]
#define POTFIT_THRMO_98 POTFIT_THRMO_97,names[97],errs[97]
#define POTFIT_THRMO_99 POTFIT_THRMO_98,names[98],errs[98]
#define POTFIT_THRMO_100 POTFIT_THRMO_99,names[99],errs[99]
#define POTFIT_THRMO_101 POTFIT_THRMO_100,names[100],errs[100]
#define POTFIT_THRMO_102 POTFIT_THRMO_101,names[101],errs[101]
#define POTFIT_THRMO_103 POTFIT_THRMO_102,names[102],errs[102]
#define POTFIT_THRMO_104 POTFIT_THRMO_103,names[103],errs[103]
#define POTFIT_THRMO_105 POTFIT_THRMO_104,names[104],errs[104]
#define POTFIT_THRMO_106 POTFIT_THRMO_105,names[105],errs[105]
#define POTFIT_THRMO_107 POTFIT_THRMO_106,names[106],errs[106]
#define POTFIT_THRMO_108 POTFIT_THRMO_107,names[107],errs[107]
#define POTFIT_THRMO_109 POTFIT_THRMO_108,names[108],errs[108]
#define POTFIT_THRMO_110 POTFIT_THRMO_109,names[109],errs[109]
#define POTFIT_THRMO_111 POTFIT_THRMO_110,names[110],errs[110]
#define POTFIT_THRMO_112 POTFIT_THRMO_111,names[111],errs[111]
#define POTFIT_THRMO_113 POTFIT_THRMO_112,names[112],errs[112]
#define POTFIT_THRMO_114 POTFIT_THRMO_113,names[113],errs[113]
#define POTFIT_THRMO_115 POTFIT_THRMO_114,names[114],errs[114]
#define POTFIT_THRMO_116 POTFIT_THRMO_115,names[115],errs[115]
#define POTFIT_THRMO_117 POTFIT_THRMO_116,names[116],errs[116]
#define POTFIT_THRMO_118 POTFIT_THRMO_117,names[117],errs[117]
#define POTFIT_THRMO_119 POTFIT_THRMO_118,names[118],errs[118]
#define POTFIT_THRMO_120 POTFIT_THRMO_119,names[119],errs[119]
#define POTFIT_THRMO_121 POTFIT_THRMO_120,names[120],errs[120]
#define POTFIT_THRMO_122 POTFIT_THRMO_121,names[121],errs[121]
#define POTFIT_THRMO_123 POTFIT_THRMO_122,names[122],errs[122]
#define POTFIT_THRMO_124 POTFIT_THRMO_123,names[123],errs[123]
#define POTFIT_THRMO_125 POTFIT_THRMO_124,names[124],errs[124]
#define POTFIT_THRMO_126 POTFIT_THRMO_125,names[125],errs[125]
#define POTFIT_THRMO_127 POTFIT_THRMO_126,names[126],errs[126]
#define POTFIT_THRMO_128 POTFIT_THRMO_127,names[127],errs[127]
#define POTFIT_THRMO_129 POTFIT_THRMO_128,names[128],errs[128]
#define POTFIT_THRMO_130 POTFIT_THRMO_129,names[129],errs[129]
#define POTFIT_THRMO_131 POTFIT_THRMO_130,names[130],errs[130]
#define POTFIT_THRMO_132 POTFIT_THRMO_131,names[131],errs[131]
#define POTFIT_THRMO_133 POTFIT_THRMO_132,names[132],errs[132]
#define POTFIT_THRMO_134 POTFIT_THRMO_133,names[133],errs[133]
#define POTFIT_THRMO_135 POTFIT_THRMO_134,names[134],errs[134]
#define POTFIT_THRMO_136 POTFIT_THRMO_135,names[135],errs[135]
#define POTFIT_THRMO_137 POTFIT_THRMO_136,names[136],errs[136]
#define POTFIT_THRMO_138 POTFIT_THRMO_137,names[137],errs[137]
#define POTFIT_THRMO_139 POTFIT_THRMO_138,names[138],errs[138]
#define POTFIT_THRMO_140 POTFIT_THRMO_139,names[139],errs[139]
#define POTFIT_THRMO_141 POTFIT_THRMO_140,names[140],errs[140]
#define POTFIT_THRMO_142 POTFIT_THRMO_141,names[141],errs[141]
#define POTFIT_THRMO_143 POTFIT_THRMO_142,names[142],errs[142]
#define POTFIT_THRMO_144 POTFIT_THRMO_143,names[143],errs[143]
#define POTFIT_THRMO_145 POTFIT_THRMO_144,names[144],errs[144]
#define POTFIT_THRMO_146 POTFIT_THRMO_145,names[145],errs[145]
#define POTFIT_THRMO_147 POTFIT_THRMO_146,names[146],errs[146]
#define POTFIT_THRMO_148 POTFIT_THRMO_147,names[147],errs[147]
#define POTFIT_THRMO_149 POTFIT_THRMO_148,names[148],errs[148]
#define POTFIT_THRMO_150 POTFIT_THRMO_149,names[149],errs[149]
#define POTFIT_THRMO_151 POTFIT_THRMO_150,names[150],errs[150]
#define POTFIT_THRMO_152 POTFIT_THRMO_151,names[151],errs[151]
#define POTFIT_THRMO_153 POTFIT_THRMO_152,names[152],errs[152]
#define POTFIT_THRMO_154 POTFIT_THRMO_153,names[153],errs[153]
#define POTFIT_THRMO_155 POTFIT_THRMO_154,names[154],errs[154]
#define POTFIT_THRMO_156 POTFIT_THRMO_155,names[155],errs[155]
#define POTFIT_THRMO_157 POTFIT_THRMO_156,names[156],errs[156]
#define POTFIT_THRMO_158 POTFIT_THRMO_157,names[157],errs[157]
#define POTFIT_THRMO_159 POTFIT_THRMO_158,names[158],errs[158]
#define POTFIT_THRMO_160 POTFIT_THRMO_159,names[159],errs[159]
#define POTFIT_THRMO_161 POTFIT_THRMO_160,names[160],errs[160]
#define POTFIT_THRMO_162 POTFIT_THRMO_161,names[161],errs[161]
#define POTFIT_THRMO_163 POTFIT_THRMO_162,names[162],errs[162]
#define POTFIT_THRMO_164 POTFIT_THRMO_163,names[163],errs[163]
#define POTFIT_THRMO_165 POTFIT_THRMO_164,names[164],errs[164]
#define POTFIT_THRMO_166 POTFIT_THRMO_165,names[165],errs[165]
#define POTFIT_THRMO_167 POTFIT_THRMO_166,names[166],errs[166]
#define POTFIT_THRMO_168 POTFIT_THRMO_167,names[167],errs[167]
#define POTFIT_THRMO_169 POTFIT_THRMO_168,names[168],errs[168]
#define POTFIT_THRMO_170 POTFIT_THRMO_169,names[169],errs[169]
#define POTFIT_THRMO_171 POTFIT_THRMO_170,names[170],errs[170]
#define POTFIT_THRMO_172 POTFIT_THRMO_171,names[171],errs[171]
#define POTFIT_THRMO_173 POTFIT_THRMO_172,names[172],errs[172]
#define POTFIT_THRMO_174 POTFIT_THRMO_173,names[173],errs[173]
#define POTFIT_THRMO_175 POTFIT_THRMO_174,names[174],errs[174]
#define POTFIT_THRMO_176 POTFIT_THRMO_175,names[175],errs[175]
#define POTFIT_THRMO_177 POTFIT_THRMO_176,names[176],errs[176]
#define POTFIT_THRMO_178 POTFIT_THRMO_177,names[177],errs[177]
#define POTFIT_THRMO_179 POTFIT_THRMO_178,names[178],errs[178]
#define POTFIT_THRMO_180 POTFIT_THRMO_179,names[179],errs[179]
#define POTFIT_THRMO_181 POTFIT_THRMO_180,names[180],errs[180]
#define POTFIT_THRMO_182 POTFIT_THRMO_181,names[181],errs[181]
#define POTFIT_THRMO_183 POTFIT_THRMO_182,names[182],errs[182]
#define POTFIT_THRMO_184 POTFIT_THRMO_183,names[183],errs[183]
#define POTFIT_THRMO_185 POTFIT_THRMO_184,names[184],errs[184]
#define POTFIT_THRMO_186 POTFIT_THRMO_185,names[185],errs[185]
#define POTFIT_THRMO_187 POTFIT_THRMO_186,names[186],errs[186]
#define POTFIT_THRMO_188 POTFIT_THRMO_187,names[187],errs[187]
#define POTFIT_THRMO_189 POTFIT_THRMO_188,names[188],errs[188]
#define POTFIT_THRMO_190 POTFIT_THRMO_189,names[189],errs[189]
#define POTFIT_THRMO_191 POTFIT_THRMO_190,names[190],errs[190]
#define POTFIT_THRMO_192 POTFIT_THRMO_191,names[191],errs[191]
#define POTFIT_THRMO_193 POTFIT_THRMO_192,names[192],errs[192]
#define POTFIT_THRMO_194 POTFIT_THRMO_193,names[193],errs[193]
#define POTFIT_THRMO_195 POTFIT_THRMO_194,names[194],errs[194]
#define POTFIT_THRMO_196 POTFIT_THRMO_195,names[195],errs[195]
#define POTFIT_THRMO_197 POTFIT_THRMO_196,names[196],errs[196]
#define POTFIT_THRMO_198 POTFIT_THRMO_197,names[197],errs[197]
#define POTFIT_THRMO_199 POTFIT_THRMO_198,names[198],errs[198]
#define POTFIT_THRMO_200 POTFIT_THRMO_199,names[199],errs[199]
#define POTFIT_THRMO_201 POTFIT_THRMO_200,names[200],errs[200]
#define POTFIT_THRMO_202 POTFIT_THRMO_201,names[201],errs[201]
#define POTFIT_THRMO_203 POTFIT_THRMO_202,names[202],errs[202]
#define POTFIT_THRMO_204 POTFIT_THRMO_203,names[203],errs[203]
#define POTFIT_THRMO_205 POTFIT_THRMO_204,names[204],errs[204]
#define POTFIT_THRMO_206 POTFIT_THRMO_205,names[205],errs[205]
#define POTFIT_THRMO_207 POTFIT_THRMO_206,names[206],errs[206]
#define POTFIT_THRMO_208 POTFIT_THRMO_207,names[207],errs[207]
#define POTFIT_THRMO_209 POTFIT_THRMO_208,names[208],errs[208]
#define POTFIT_THRMO_210 POTFIT_THRMO_209,names[209],errs[209]
#define POTFIT_THRMO_211 POTFIT_THRMO_210,names[210],errs[210]
#define POTFIT_THRMO_212 POTFIT_THRMO_211,names[211],errs[211]
#define POTFIT_THRMO_213 POTFIT_THRMO_212,names[212],errs[212]
#define POTFIT_THRMO_214 POTFIT_THRMO_213,names[213],errs[213]
#define POTFIT_THRMO_215 POTFIT_THRMO_214,names[214],errs[214]
#define POTFIT_THRMO_216 POTFIT_THRMO_215,names[215],errs[215]
#define POTFIT_THRMO_217 POTFIT_THRMO_216,names[216],errs[216]
#define POTFIT_THRMO_218 POTFIT_THRMO_217,names[217],errs[217]
#define POTFIT_THRMO_219 POTFIT_THRMO_218,names[218],errs[218]
#define POTFIT_THRMO_220 POTFIT_THRMO_219,names[219],errs[219]
#define POTFIT_THRMO_221 POTFIT_THRMO_220,names[220],errs[220]
#define POTFIT_THRMO_222 POTFIT_THRMO_221,names[221],errs[221]
#define POTFIT_THRMO_223 POTFIT_THRMO_222,names[222],errs[222]
#define POTFIT_THRMO_224 POTFIT_THRMO_223,names[223],errs[223]
#define POTFIT_THRMO_225 POTFIT_THRMO_224,names[224],errs[224]
#define POTFIT_THRMO_226 POTFIT_THRMO_225,names[225],errs[225]
#define POTFIT_THRMO_227 POTFIT_THRMO_226,names[226],errs[226]
#define POTFIT_THRMO_228 POTFIT_THRMO_227,names[227],errs[227]
#define POTFIT_THRMO_229 POTFIT_THRMO_228,names[228],errs[228]
#define POTFIT_THRMO_230 POTFIT_THRMO_229,names[229],errs[229]
#define POTFIT_THRMO_231 POTFIT_THRMO_230,names[230],errs[230]
#define POTFIT_THRMO_232 POTFIT_THRMO_231,names[231],errs[231]
#define POTFIT_THRMO_233 POTFIT_THRMO_232,names[232],errs[232]
#define POTFIT_THRMO_234 POTFIT_THRMO_233,names[233],errs[233]
#define POTFIT_THRMO_235 POTFIT_THRMO_234,names[234],errs[234]
#define POTFIT_THRMO_236 POTFIT_THRMO_235,names[235],errs[235]
#define POTFIT_THRMO_237 POTFIT_THRMO_236,names[236],errs[236]
#define POTFIT_THRMO_238 POTFIT_THRMO_237,names[237],errs[237]
#define POTFIT_THRMO_239 POTFIT_THRMO_238,names[238],errs[238]
#define POTFIT_THRMO_240 POTFIT_THRMO_239,names[239],errs[239]
#define POTFIT_THRMO_241 POTFIT_THRMO_240,names[240],errs[240]
#define POTFIT_THRMO_242 POTFIT_THRMO_241,names[241],errs[241]
#define POTFIT_THRMO_243 POTFIT_THRMO_242,names[242],errs[242]
#define POTFIT_THRMO_244 POTFIT_THRMO_243,names[243],errs[243]
#define POTFIT_THRMO_245 POTFIT_THRMO_244,names[244],errs[244]
#define POTFIT_THRMO_246 POTFIT_THRMO_245,names[245],errs[245]
#define POTFIT_THRMO_247 POTFIT_THRMO_246,names[246],errs[246]
#define POTFIT_THRMO_248 POTFIT_THRMO_247,names[247],errs[247]
#define POTFIT_THRMO_249 POTFIT_THRMO_248,names[248],errs[248]
#define POTFIT_THRMO_250 POTFIT_THRMO_249,names[249],errs[249]
#define POTFIT_THRMO_251 POTFIT_THRMO_250,names[250],errs[250]
#define POTFIT_THRMO_252 POTFIT_THRMO_251,names[251],errs[251]
#define POTFIT_THRMO_253 POTFIT_THRMO_252,names[252],errs[252]
#define POTFIT_THRMO_254 POTFIT_THRMO_253,names[253],errs[253]
#define POTFIT_THRMO_255 POTFIT_THRMO_254,names[254],errs[254]
#define POTFIT_THRMO_256 POTFIT_THRMO_255,names[255],errs[255]
#define RET_POTFIT_THRMO_0 if(nconfigs==0) return ThermoDynamics(POTFIT_THRMO_0); RET_POTFIT_THRMO_1
#define RET_POTFIT_THRMO_1 if(nconfigs==1) return ThermoDynamics(POTFIT_THRMO_1); RET_POTFIT_THRMO_2
#define RET_POTFIT_THRMO_2 if(nconfigs==2) return ThermoDynamics(POTFIT_THRMO_2); RET_POTFIT_THRMO_3
#define RET_POTFIT_THRMO_3 if(nconfigs==3) return ThermoDynamics(POTFIT_THRMO_3); RET_POTFIT_THRMO_4
#define RET_POTFIT_THRMO_4 if(nconfigs==4) return ThermoDynamics(POTFIT_THRMO_4); RET_POTFIT_THRMO_5
#define RET_POTFIT_THRMO_5 if(nconfigs==5) return ThermoDynamics(POTFIT_THRMO_5); RET_POTFIT_THRMO_6
#define RET_POTFIT_THRMO_6 if(nconfigs==6) return ThermoDynamics(POTFIT_THRMO_6); RET_POTFIT_THRMO_7
#define RET_POTFIT_THRMO_7 if(nconfigs==7) return ThermoDynamics(POTFIT_THRMO_7); RET_POTFIT_THRMO_8
#define RET_POTFIT_THRMO_8 if(nconfigs==8) return ThermoDynamics(POTFIT_THRMO_8); RET_POTFIT_THRMO_9
#define RET_POTFIT_THRMO_9 if(nconfigs==9) return ThermoDynamics(POTFIT_THRMO_9); RET_POTFIT_THRMO_10
#define RET_POTFIT_THRMO_10 if(nconfigs==10) return ThermoDynamics(POTFIT_THRMO_10); RET_POTFIT_THRMO_11
#define RET_POTFIT_THRMO_11 if(nconfigs==11) return ThermoDynamics(POTFIT_THRMO_11); RET_POTFIT_THRMO_12
#define RET_POTFIT_THRMO_12 if(nconfigs==12) return ThermoDynamics(POTFIT_THRMO_12); RET_POTFIT_THRMO_13
#define RET_POTFIT_THRMO_13 if(nconfigs==13) return ThermoDynamics(POTFIT_THRMO_13); RET_POTFIT_THRMO_14
#define RET_POTFIT_THRMO_14 if(nconfigs==14) return ThermoDynamics(POTFIT_THRMO_14); RET_POTFIT_THRMO_15
#define RET_POTFIT_THRMO_15 if(nconfigs==15) return ThermoDynamics(POTFIT_THRMO_15); RET_POTFIT_THRMO_16
#define RET_POTFIT_THRMO_16 if(nconfigs==16) return ThermoDynamics(POTFIT_THRMO_16); RET_POTFIT_THRMO_17
#define RET_POTFIT_THRMO_17 if(nconfigs==17) return ThermoDynamics(POTFIT_THRMO_17); RET_POTFIT_THRMO_18
#define RET_POTFIT_THRMO_18 if(nconfigs==18) return ThermoDynamics(POTFIT_THRMO_18); RET_POTFIT_THRMO_19
#define RET_POTFIT_THRMO_19 if(nconfigs==19) return ThermoDynamics(POTFIT_THRMO_19); RET_POTFIT_THRMO_20
#define RET_POTFIT_THRMO_20 if(nconfigs==20) return ThermoDynamics(POTFIT_THRMO_20); RET_POTFIT_THRMO_21
#define RET_POTFIT_THRMO_21 if(nconfigs==21) return ThermoDynamics(POTFIT_THRMO_21); RET_POTFIT_THRMO_22
#define RET_POTFIT_THRMO_22 if(nconfigs==22) return ThermoDynamics(POTFIT_THRMO_22); RET_POTFIT_THRMO_23
#define RET_POTFIT_THRMO_23 if(nconfigs==23) return ThermoDynamics(POTFIT_THRMO_23); RET_POTFIT_THRMO_24
#define RET_POTFIT_THRMO_24 if(nconfigs==24) return ThermoDynamics(POTFIT_THRMO_24); RET_POTFIT_THRMO_25
#define RET_POTFIT_THRMO_25 if(nconfigs==25) return ThermoDynamics(POTFIT_THRMO_25); RET_POTFIT_THRMO_26
#define RET_POTFIT_THRMO_26 if(nconfigs==26) return ThermoDynamics(POTFIT_THRMO_26); RET_POTFIT_THRMO_27
#define RET_POTFIT_THRMO_27 if(nconfigs==27) return ThermoDynamics(POTFIT_THRMO_27); RET_POTFIT_THRMO_28
#define RET_POTFIT_THRMO_28 if(nconfigs==28) return ThermoDynamics(POTFIT_THRMO_28); RET_POTFIT_THRMO_29
#define RET_POTFIT_THRMO_29 if(nconfigs==29) return ThermoDynamics(POTFIT_THRMO_29); RET_POTFIT_THRMO_30
#define RET_POTFIT_THRMO_30 if(nconfigs==30) return ThermoDynamics(POTFIT_THRMO_30); RET_POTFIT_THRMO_31
#define RET_POTFIT_THRMO_31 if(nconfigs==31) return ThermoDynamics(POTFIT_THRMO_31); RET_POTFIT_THRMO_32
#define RET_POTFIT_THRMO_32 if(nconfigs==32) return ThermoDynamics(POTFIT_THRMO_32); RET_POTFIT_THRMO_33
#define RET_POTFIT_THRMO_33 if(nconfigs==33) return ThermoDynamics(POTFIT_THRMO_33); RET_POTFIT_THRMO_34
#define RET_POTFIT_THRMO_34 if(nconfigs==34) return ThermoDynamics(POTFIT_THRMO_34); RET_POTFIT_THRMO_35
#define RET_POTFIT_THRMO_35 if(nconfigs==35) return ThermoDynamics(POTFIT_THRMO_35); RET_POTFIT_THRMO_36
#define RET_POTFIT_THRMO_36 if(nconfigs==36) return ThermoDynamics(POTFIT_THRMO_36); RET_POTFIT_THRMO_37
#define RET_POTFIT_THRMO_37 if(nconfigs==37) return ThermoDynamics(POTFIT_THRMO_37); RET_POTFIT_THRMO_38
#define RET_POTFIT_THRMO_38 if(nconfigs==38) return ThermoDynamics(POTFIT_THRMO_38); RET_POTFIT_THRMO_39
#define RET_POTFIT_THRMO_39 if(nconfigs==39) return ThermoDynamics(POTFIT_THRMO_39); RET_POTFIT_THRMO_40
#define RET_POTFIT_THRMO_40 if(nconfigs==40) return ThermoDynamics(POTFIT_THRMO_40); RET_POTFIT_THRMO_41
#define RET_POTFIT_THRMO_41 if(nconfigs==41) return ThermoDynamics(POTFIT_THRMO_41); RET_POTFIT_THRMO_42
#define RET_POTFIT_THRMO_42 if(nconfigs==42) return ThermoDynamics(POTFIT_THRMO_42); RET_POTFIT_THRMO_43
#define RET_POTFIT_THRMO_43 if(nconfigs==43) return ThermoDynamics(POTFIT_THRMO_43); RET_POTFIT_THRMO_44
#define RET_POTFIT_THRMO_44 if(nconfigs==44) return ThermoDynamics(POTFIT_THRMO_44); RET_POTFIT_THRMO_45
#define RET_POTFIT_THRMO_45 if(nconfigs==45) return ThermoDynamics(POTFIT_THRMO_45); RET_POTFIT_THRMO_46
#define RET_POTFIT_THRMO_46 if(nconfigs==46) return ThermoDynamics(POTFIT_THRMO_46); RET_POTFIT_THRMO_47
#define RET_POTFIT_THRMO_47 if(nconfigs==47) return ThermoDynamics(POTFIT_THRMO_47); RET_POTFIT_THRMO_48
#define RET_POTFIT_THRMO_48 if(nconfigs==48) return ThermoDynamics(POTFIT_THRMO_48); RET_POTFIT_THRMO_49
#define RET_POTFIT_THRMO_49 if(nconfigs==49) return ThermoDynamics(POTFIT_THRMO_49); RET_POTFIT_THRMO_50
#define RET_POTFIT_THRMO_50 if(nconfigs==50) return ThermoDynamics(POTFIT_THRMO_50); RET_POTFIT_THRMO_51
#define RET_POTFIT_THRMO_51 if(nconfigs==51) return ThermoDynamics(POTFIT_THRMO_51); RET_POTFIT_THRMO_52
#define RET_POTFIT_THRMO_52 if(nconfigs==52) return ThermoDynamics(POTFIT_THRMO_52); RET_POTFIT_THRMO_53
#define RET_POTFIT_THRMO_53 if(nconfigs==53) return ThermoDynamics(POTFIT_THRMO_53); RET_POTFIT_THRMO_54
#define RET_POTFIT_THRMO_54 if(nconfigs==54) return ThermoDynamics(POTFIT_THRMO_54); RET_POTFIT_THRMO_55
#define RET_POTFIT_THRMO_55 if(nconfigs==55) return ThermoDynamics(POTFIT_THRMO_55); RET_POTFIT_THRMO_56
#define RET_POTFIT_THRMO_56 if(nconfigs==56) return ThermoDynamics(POTFIT_THRMO_56); RET_POTFIT_THRMO_57
#define RET_POTFIT_THRMO_57 if(nconfigs==57) return ThermoDynamics(POTFIT_THRMO_57); RET_POTFIT_THRMO_58
#define RET_POTFIT_THRMO_58 if(nconfigs==58) return ThermoDynamics(POTFIT_THRMO_58); RET_POTFIT_THRMO_59
#define RET_POTFIT_THRMO_59 if(nconfigs==59) return ThermoDynamics(POTFIT_THRMO_59); RET_POTFIT_THRMO_60
#define RET_POTFIT_THRMO_60 if(nconfigs==60) return ThermoDynamics(POTFIT_THRMO_60); RET_POTFIT_THRMO_61
#define RET_POTFIT_THRMO_61 if(nconfigs==61) return ThermoDynamics(POTFIT_THRMO_61); RET_POTFIT_THRMO_62
#define RET_POTFIT_THRMO_62 if(nconfigs==62) return ThermoDynamics(POTFIT_THRMO_62); RET_POTFIT_THRMO_63
#define RET_POTFIT_THRMO_63 if(nconfigs==63) return ThermoDynamics(POTFIT_THRMO_63); RET_POTFIT_THRMO_64
#define RET_POTFIT_THRMO_64 if(nconfigs==64) return ThermoDynamics(POTFIT_THRMO_64); RET_POTFIT_THRMO_65
#define RET_POTFIT_THRMO_65 if(nconfigs==65) return ThermoDynamics(POTFIT_THRMO_65); RET_POTFIT_THRMO_66
#define RET_POTFIT_THRMO_66 if(nconfigs==66) return ThermoDynamics(POTFIT_THRMO_66); RET_POTFIT_THRMO_67
#define RET_POTFIT_THRMO_67 if(nconfigs==67) return ThermoDynamics(POTFIT_THRMO_67); RET_POTFIT_THRMO_68
#define RET_POTFIT_THRMO_68 if(nconfigs==68) return ThermoDynamics(POTFIT_THRMO_68); RET_POTFIT_THRMO_69
#define RET_POTFIT_THRMO_69 if(nconfigs==69) return ThermoDynamics(POTFIT_THRMO_69); RET_POTFIT_THRMO_70
#define RET_POTFIT_THRMO_70 if(nconfigs==70) return ThermoDynamics(POTFIT_THRMO_70); RET_POTFIT_THRMO_71
#define RET_POTFIT_THRMO_71 if(nconfigs==71) return ThermoDynamics(POTFIT_THRMO_71); RET_POTFIT_THRMO_72
#define RET_POTFIT_THRMO_72 if(nconfigs==72) return ThermoDynamics(POTFIT_THRMO_72); RET_POTFIT_THRMO_73
#define RET_POTFIT_THRMO_73 if(nconfigs==73) return ThermoDynamics(POTFIT_THRMO_73); RET_POTFIT_THRMO_74
#define RET_POTFIT_THRMO_74 if(nconfigs==74) return ThermoDynamics(POTFIT_THRMO_74); RET_POTFIT_THRMO_75
#define RET_POTFIT_THRMO_75 if(nconfigs==75) return ThermoDynamics(POTFIT_THRMO_75); RET_POTFIT_THRMO_76
#define RET_POTFIT_THRMO_76 if(nconfigs==76) return ThermoDynamics(POTFIT_THRMO_76); RET_POTFIT_THRMO_77
#define RET_POTFIT_THRMO_77 if(nconfigs==77) return ThermoDynamics(POTFIT_THRMO_77); RET_POTFIT_THRMO_78
#define RET_POTFIT_THRMO_78 if(nconfigs==78) return ThermoDynamics(POTFIT_THRMO_78); RET_POTFIT_THRMO_79
#define RET_POTFIT_THRMO_79 if(nconfigs==79) return ThermoDynamics(POTFIT_THRMO_79); RET_POTFIT_THRMO_80
#define RET_POTFIT_THRMO_80 if(nconfigs==80) return ThermoDynamics(POTFIT_THRMO_80); RET_POTFIT_THRMO_81
#define RET_POTFIT_THRMO_81 if(nconfigs==81) return ThermoDynamics(POTFIT_THRMO_81); RET_POTFIT_THRMO_82
#define RET_POTFIT_THRMO_82 if(nconfigs==82) return ThermoDynamics(POTFIT_THRMO_82); RET_POTFIT_THRMO_83
#define RET_POTFIT_THRMO_83 if(nconfigs==83) return ThermoDynamics(POTFIT_THRMO_83); RET_POTFIT_THRMO_84
#define RET_POTFIT_THRMO_84 if(nconfigs==84) return ThermoDynamics(POTFIT_THRMO_84); RET_POTFIT_THRMO_85
#define RET_POTFIT_THRMO_85 if(nconfigs==85) return ThermoDynamics(POTFIT_THRMO_85); RET_POTFIT_THRMO_86
#define RET_POTFIT_THRMO_86 if(nconfigs==86) return ThermoDynamics(POTFIT_THRMO_86); RET_POTFIT_THRMO_87
#define RET_POTFIT_THRMO_87 if(nconfigs==87) return ThermoDynamics(POTFIT_THRMO_87); RET_POTFIT_THRMO_88
#define RET_POTFIT_THRMO_88 if(nconfigs==88) return ThermoDynamics(POTFIT_THRMO_88); RET_POTFIT_THRMO_89
#define RET_POTFIT_THRMO_89 if(nconfigs==89) return ThermoDynamics(POTFIT_THRMO_89); RET_POTFIT_THRMO_90
#define RET_POTFIT_THRMO_90 if(nconfigs==90) return ThermoDynamics(POTFIT_THRMO_90); RET_POTFIT_THRMO_91
#define RET_POTFIT_THRMO_91 if(nconfigs==91) return ThermoDynamics(POTFIT_THRMO_91); RET_POTFIT_THRMO_92
#define RET_POTFIT_THRMO_92 if(nconfigs==92) return ThermoDynamics(POTFIT_THRMO_92); RET_POTFIT_THRMO_93
#define RET_POTFIT_THRMO_93 if(nconfigs==93) return ThermoDynamics(POTFIT_THRMO_93); RET_POTFIT_THRMO_94
#define RET_POTFIT_THRMO_94 if(nconfigs==94) return ThermoDynamics(POTFIT_THRMO_94); RET_POTFIT_THRMO_95
#define RET_POTFIT_THRMO_95 if(nconfigs==95) return ThermoDynamics(POTFIT_THRMO_95); RET_POTFIT_THRMO_96
#define RET_POTFIT_THRMO_96 if(nconfigs==96) return ThermoDynamics(POTFIT_THRMO_96); RET_POTFIT_THRMO_97
#define RET_POTFIT_THRMO_97 if(nconfigs==97) return ThermoDynamics(POTFIT_THRMO_97); RET_POTFIT_THRMO_98
#define RET_POTFIT_THRMO_98 if(nconfigs==98) return ThermoDynamics(POTFIT_THRMO_98); RET_POTFIT_THRMO_99
#define RET_POTFIT_THRMO_99 if(nconfigs==99) return ThermoDynamics(POTFIT_THRMO_99); RET_POTFIT_THRMO_100
#define RET_POTFIT_THRMO_100 if(nconfigs==100) return ThermoDynamics(POTFIT_THRMO_100); RET_POTFIT_THRMO_101
#define RET_POTFIT_THRMO_101 if(nconfigs==101) return ThermoDynamics(POTFIT_THRMO_101); RET_POTFIT_THRMO_102
#define RET_POTFIT_THRMO_102 if(nconfigs==102) return ThermoDynamics(POTFIT_THRMO_102); RET_POTFIT_THRMO_103
#define RET_POTFIT_THRMO_103 if(nconfigs==103) return ThermoDynamics(POTFIT_THRMO_103); RET_POTFIT_THRMO_104
#define RET_POTFIT_THRMO_104 if(nconfigs==104) return ThermoDynamics(POTFIT_THRMO_104); RET_POTFIT_THRMO_105
#define RET_POTFIT_THRMO_105 if(nconfigs==105) return ThermoDynamics(POTFIT_THRMO_105); RET_POTFIT_THRMO_106
#define RET_POTFIT_THRMO_106 if(nconfigs==106) return ThermoDynamics(POTFIT_THRMO_106); RET_POTFIT_THRMO_107
#define RET_POTFIT_THRMO_107 if(nconfigs==107) return ThermoDynamics(POTFIT_THRMO_107); RET_POTFIT_THRMO_108
#define RET_POTFIT_THRMO_108 if(nconfigs==108) return ThermoDynamics(POTFIT_THRMO_108); RET_POTFIT_THRMO_109
#define RET_POTFIT_THRMO_109 if(nconfigs==109) return ThermoDynamics(POTFIT_THRMO_109); RET_POTFIT_THRMO_110
#define RET_POTFIT_THRMO_110 if(nconfigs==110) return ThermoDynamics(POTFIT_THRMO_110); RET_POTFIT_THRMO_111
#define RET_POTFIT_THRMO_111 if(nconfigs==111) return ThermoDynamics(POTFIT_THRMO_111); RET_POTFIT_THRMO_112
#define RET_POTFIT_THRMO_112 if(nconfigs==112) return ThermoDynamics(POTFIT_THRMO_112); RET_POTFIT_THRMO_113
#define RET_POTFIT_THRMO_113 if(nconfigs==113) return ThermoDynamics(POTFIT_THRMO_113); RET_POTFIT_THRMO_114
#define RET_POTFIT_THRMO_114 if(nconfigs==114) return ThermoDynamics(POTFIT_THRMO_114); RET_POTFIT_THRMO_115
#define RET_POTFIT_THRMO_115 if(nconfigs==115) return ThermoDynamics(POTFIT_THRMO_115); RET_POTFIT_THRMO_116
#define RET_POTFIT_THRMO_116 if(nconfigs==116) return ThermoDynamics(POTFIT_THRMO_116); RET_POTFIT_THRMO_117
#define RET_POTFIT_THRMO_117 if(nconfigs==117) return ThermoDynamics(POTFIT_THRMO_117); RET_POTFIT_THRMO_118
#define RET_POTFIT_THRMO_118 if(nconfigs==118) return ThermoDynamics(POTFIT_THRMO_118); RET_POTFIT_THRMO_119
#define RET_POTFIT_THRMO_119 if(nconfigs==119) return ThermoDynamics(POTFIT_THRMO_119); RET_POTFIT_THRMO_120
#define RET_POTFIT_THRMO_120 if(nconfigs==120) return ThermoDynamics(POTFIT_THRMO_120); RET_POTFIT_THRMO_121
#define RET_POTFIT_THRMO_121 if(nconfigs==121) return ThermoDynamics(POTFIT_THRMO_121); RET_POTFIT_THRMO_122
#define RET_POTFIT_THRMO_122 if(nconfigs==122) return ThermoDynamics(POTFIT_THRMO_122); RET_POTFIT_THRMO_123
#define RET_POTFIT_THRMO_123 if(nconfigs==123) return ThermoDynamics(POTFIT_THRMO_123); RET_POTFIT_THRMO_124
#define RET_POTFIT_THRMO_124 if(nconfigs==124) return ThermoDynamics(POTFIT_THRMO_124); RET_POTFIT_THRMO_125
#define RET_POTFIT_THRMO_125 if(nconfigs==125) return ThermoDynamics(POTFIT_THRMO_125); RET_POTFIT_THRMO_126
#define RET_POTFIT_THRMO_126 if(nconfigs==126) return ThermoDynamics(POTFIT_THRMO_126); RET_POTFIT_THRMO_127
#define RET_POTFIT_THRMO_127 if(nconfigs==127) return ThermoDynamics(POTFIT_THRMO_127); RET_POTFIT_THRMO_128
#define RET_POTFIT_THRMO_128 if(nconfigs==128) return ThermoDynamics(POTFIT_THRMO_128); RET_POTFIT_THRMO_129
#define RET_POTFIT_THRMO_129 if(nconfigs==129) return ThermoDynamics(POTFIT_THRMO_129); RET_POTFIT_THRMO_130
#define RET_POTFIT_THRMO_130 if(nconfigs==130) return ThermoDynamics(POTFIT_THRMO_130); RET_POTFIT_THRMO_131
#define RET_POTFIT_THRMO_131 if(nconfigs==131) return ThermoDynamics(POTFIT_THRMO_131); RET_POTFIT_THRMO_132
#define RET_POTFIT_THRMO_132 if(nconfigs==132) return ThermoDynamics(POTFIT_THRMO_132); RET_POTFIT_THRMO_133
#define RET_POTFIT_THRMO_133 if(nconfigs==133) return ThermoDynamics(POTFIT_THRMO_133); RET_POTFIT_THRMO_134
#define RET_POTFIT_THRMO_134 if(nconfigs==134) return ThermoDynamics(POTFIT_THRMO_134); RET_POTFIT_THRMO_135
#define RET_POTFIT_THRMO_135 if(nconfigs==135) return ThermoDynamics(POTFIT_THRMO_135); RET_POTFIT_THRMO_136
#define RET_POTFIT_THRMO_136 if(nconfigs==136) return ThermoDynamics(POTFIT_THRMO_136); RET_POTFIT_THRMO_137
#define RET_POTFIT_THRMO_137 if(nconfigs==137) return ThermoDynamics(POTFIT_THRMO_137); RET_POTFIT_THRMO_138
#define RET_POTFIT_THRMO_138 if(nconfigs==138) return ThermoDynamics(POTFIT_THRMO_138); RET_POTFIT_THRMO_139
#define RET_POTFIT_THRMO_139 if(nconfigs==139) return ThermoDynamics(POTFIT_THRMO_139); RET_POTFIT_THRMO_140
#define RET_POTFIT_THRMO_140 if(nconfigs==140) return ThermoDynamics(POTFIT_THRMO_140); RET_POTFIT_THRMO_141
#define RET_POTFIT_THRMO_141 if(nconfigs==141) return ThermoDynamics(POTFIT_THRMO_141); RET_POTFIT_THRMO_142
#define RET_POTFIT_THRMO_142 if(nconfigs==142) return ThermoDynamics(POTFIT_THRMO_142); RET_POTFIT_THRMO_143
#define RET_POTFIT_THRMO_143 if(nconfigs==143) return ThermoDynamics(POTFIT_THRMO_143); RET_POTFIT_THRMO_144
#define RET_POTFIT_THRMO_144 if(nconfigs==144) return ThermoDynamics(POTFIT_THRMO_144); RET_POTFIT_THRMO_145
#define RET_POTFIT_THRMO_145 if(nconfigs==145) return ThermoDynamics(POTFIT_THRMO_145); RET_POTFIT_THRMO_146
#define RET_POTFIT_THRMO_146 if(nconfigs==146) return ThermoDynamics(POTFIT_THRMO_146); RET_POTFIT_THRMO_147
#define RET_POTFIT_THRMO_147 if(nconfigs==147) return ThermoDynamics(POTFIT_THRMO_147); RET_POTFIT_THRMO_148
#define RET_POTFIT_THRMO_148 if(nconfigs==148) return ThermoDynamics(POTFIT_THRMO_148); RET_POTFIT_THRMO_149
#define RET_POTFIT_THRMO_149 if(nconfigs==149) return ThermoDynamics(POTFIT_THRMO_149); RET_POTFIT_THRMO_150
#define RET_POTFIT_THRMO_150 if(nconfigs==150) return ThermoDynamics(POTFIT_THRMO_150); RET_POTFIT_THRMO_151
#define RET_POTFIT_THRMO_151 if(nconfigs==151) return ThermoDynamics(POTFIT_THRMO_151); RET_POTFIT_THRMO_152
#define RET_POTFIT_THRMO_152 if(nconfigs==152) return ThermoDynamics(POTFIT_THRMO_152); RET_POTFIT_THRMO_153
#define RET_POTFIT_THRMO_153 if(nconfigs==153) return ThermoDynamics(POTFIT_THRMO_153); RET_POTFIT_THRMO_154
#define RET_POTFIT_THRMO_154 if(nconfigs==154) return ThermoDynamics(POTFIT_THRMO_154); RET_POTFIT_THRMO_155
#define RET_POTFIT_THRMO_155 if(nconfigs==155) return ThermoDynamics(POTFIT_THRMO_155); RET_POTFIT_THRMO_156
#define RET_POTFIT_THRMO_156 if(nconfigs==156) return ThermoDynamics(POTFIT_THRMO_156); RET_POTFIT_THRMO_157
#define RET_POTFIT_THRMO_157 if(nconfigs==157) return ThermoDynamics(POTFIT_THRMO_157); RET_POTFIT_THRMO_158
#define RET_POTFIT_THRMO_158 if(nconfigs==158) return ThermoDynamics(POTFIT_THRMO_158); RET_POTFIT_THRMO_159
#define RET_POTFIT_THRMO_159 if(nconfigs==159) return ThermoDynamics(POTFIT_THRMO_159); RET_POTFIT_THRMO_160
#define RET_POTFIT_THRMO_160 if(nconfigs==160) return ThermoDynamics(POTFIT_THRMO_160); RET_POTFIT_THRMO_161
#define RET_POTFIT_THRMO_161 if(nconfigs==161) return ThermoDynamics(POTFIT_THRMO_161); RET_POTFIT_THRMO_162
#define RET_POTFIT_THRMO_162 if(nconfigs==162) return ThermoDynamics(POTFIT_THRMO_162); RET_POTFIT_THRMO_163
#define RET_POTFIT_THRMO_163 if(nconfigs==163) return ThermoDynamics(POTFIT_THRMO_163); RET_POTFIT_THRMO_164
#define RET_POTFIT_THRMO_164 if(nconfigs==164) return ThermoDynamics(POTFIT_THRMO_164); RET_POTFIT_THRMO_165
#define RET_POTFIT_THRMO_165 if(nconfigs==165) return ThermoDynamics(POTFIT_THRMO_165); RET_POTFIT_THRMO_166
#define RET_POTFIT_THRMO_166 if(nconfigs==166) return ThermoDynamics(POTFIT_THRMO_166); RET_POTFIT_THRMO_167
#define RET_POTFIT_THRMO_167 if(nconfigs==167) return ThermoDynamics(POTFIT_THRMO_167); RET_POTFIT_THRMO_168
#define RET_POTFIT_THRMO_168 if(nconfigs==168) return ThermoDynamics(POTFIT_THRMO_168); RET_POTFIT_THRMO_169
#define RET_POTFIT_THRMO_169 if(nconfigs==169) return ThermoDynamics(POTFIT_THRMO_169); RET_POTFIT_THRMO_170
#define RET_POTFIT_THRMO_170 if(nconfigs==170) return ThermoDynamics(POTFIT_THRMO_170); RET_POTFIT_THRMO_171
#define RET_POTFIT_THRMO_171 if(nconfigs==171) return ThermoDynamics(POTFIT_THRMO_171); RET_POTFIT_THRMO_172
#define RET_POTFIT_THRMO_172 if(nconfigs==172) return ThermoDynamics(POTFIT_THRMO_172); RET_POTFIT_THRMO_173
#define RET_POTFIT_THRMO_173 if(nconfigs==173) return ThermoDynamics(POTFIT_THRMO_173); RET_POTFIT_THRMO_174
#define RET_POTFIT_THRMO_174 if(nconfigs==174) return ThermoDynamics(POTFIT_THRMO_174); RET_POTFIT_THRMO_175
#define RET_POTFIT_THRMO_175 if(nconfigs==175) return ThermoDynamics(POTFIT_THRMO_175); RET_POTFIT_THRMO_176
#define RET_POTFIT_THRMO_176 if(nconfigs==176) return ThermoDynamics(POTFIT_THRMO_176); RET_POTFIT_THRMO_177
#define RET_POTFIT_THRMO_177 if(nconfigs==177) return ThermoDynamics(POTFIT_THRMO_177); RET_POTFIT_THRMO_178
#define RET_POTFIT_THRMO_178 if(nconfigs==178) return ThermoDynamics(POTFIT_THRMO_178); RET_POTFIT_THRMO_179
#define RET_POTFIT_THRMO_179 if(nconfigs==179) return ThermoDynamics(POTFIT_THRMO_179); RET_POTFIT_THRMO_180
#define RET_POTFIT_THRMO_180 if(nconfigs==180) return ThermoDynamics(POTFIT_THRMO_180); RET_POTFIT_THRMO_181
#define RET_POTFIT_THRMO_181 if(nconfigs==181) return ThermoDynamics(POTFIT_THRMO_181); RET_POTFIT_THRMO_182
#define RET_POTFIT_THRMO_182 if(nconfigs==182) return ThermoDynamics(POTFIT_THRMO_182); RET_POTFIT_THRMO_183
#define RET_POTFIT_THRMO_183 if(nconfigs==183) return ThermoDynamics(POTFIT_THRMO_183); RET_POTFIT_THRMO_184
#define RET_POTFIT_THRMO_184 if(nconfigs==184) return ThermoDynamics(POTFIT_THRMO_184); RET_POTFIT_THRMO_185
#define RET_POTFIT_THRMO_185 if(nconfigs==185) return ThermoDynamics(POTFIT_THRMO_185); RET_POTFIT_THRMO_186
#define RET_POTFIT_THRMO_186 if(nconfigs==186) return ThermoDynamics(POTFIT_THRMO_186); RET_POTFIT_THRMO_187
#define RET_POTFIT_THRMO_187 if(nconfigs==187) return ThermoDynamics(POTFIT_THRMO_187); RET_POTFIT_THRMO_188
#define RET_POTFIT_THRMO_188 if(nconfigs==188) return ThermoDynamics(POTFIT_THRMO_188); RET_POTFIT_THRMO_189
#define RET_POTFIT_THRMO_189 if(nconfigs==189) return ThermoDynamics(POTFIT_THRMO_189); RET_POTFIT_THRMO_190
#define RET_POTFIT_THRMO_190 if(nconfigs==190) return ThermoDynamics(POTFIT_THRMO_190); RET_POTFIT_THRMO_191
#define RET_POTFIT_THRMO_191 if(nconfigs==191) return ThermoDynamics(POTFIT_THRMO_191); RET_POTFIT_THRMO_192
#define RET_POTFIT_THRMO_192 if(nconfigs==192) return ThermoDynamics(POTFIT_THRMO_192); RET_POTFIT_THRMO_193
#define RET_POTFIT_THRMO_193 if(nconfigs==193) return ThermoDynamics(POTFIT_THRMO_193); RET_POTFIT_THRMO_194
#define RET_POTFIT_THRMO_194 if(nconfigs==194) return ThermoDynamics(POTFIT_THRMO_194); RET_POTFIT_THRMO_195
#define RET_POTFIT_THRMO_195 if(nconfigs==195) return ThermoDynamics(POTFIT_THRMO_195); RET_POTFIT_THRMO_196
#define RET_POTFIT_THRMO_196 if(nconfigs==196) return ThermoDynamics(POTFIT_THRMO_196); RET_POTFIT_THRMO_197
#define RET_POTFIT_THRMO_197 if(nconfigs==197) return ThermoDynamics(POTFIT_THRMO_197); RET_POTFIT_THRMO_198
#define RET_POTFIT_THRMO_198 if(nconfigs==198) return ThermoDynamics(POTFIT_THRMO_198); RET_POTFIT_THRMO_199
#define RET_POTFIT_THRMO_199 if(nconfigs==199) return ThermoDynamics(POTFIT_THRMO_199); RET_POTFIT_THRMO_200
#define RET_POTFIT_THRMO_200 if(nconfigs==200) return ThermoDynamics(POTFIT_THRMO_200); RET_POTFIT_THRMO_201
#define RET_POTFIT_THRMO_201 if(nconfigs==201) return ThermoDynamics(POTFIT_THRMO_201); RET_POTFIT_THRMO_202
#define RET_POTFIT_THRMO_202 if(nconfigs==202) return ThermoDynamics(POTFIT_THRMO_202); RET_POTFIT_THRMO_203
#define RET_POTFIT_THRMO_203 if(nconfigs==203) return ThermoDynamics(POTFIT_THRMO_203); RET_POTFIT_THRMO_204
#define RET_POTFIT_THRMO_204 if(nconfigs==204) return ThermoDynamics(POTFIT_THRMO_204); RET_POTFIT_THRMO_205
#define RET_POTFIT_THRMO_205 if(nconfigs==205) return ThermoDynamics(POTFIT_THRMO_205); RET_POTFIT_THRMO_206
#define RET_POTFIT_THRMO_206 if(nconfigs==206) return ThermoDynamics(POTFIT_THRMO_206); RET_POTFIT_THRMO_207
#define RET_POTFIT_THRMO_207 if(nconfigs==207) return ThermoDynamics(POTFIT_THRMO_207); RET_POTFIT_THRMO_208
#define RET_POTFIT_THRMO_208 if(nconfigs==208) return ThermoDynamics(POTFIT_THRMO_208); RET_POTFIT_THRMO_209
#define RET_POTFIT_THRMO_209 if(nconfigs==209) return ThermoDynamics(POTFIT_THRMO_209); RET_POTFIT_THRMO_210
#define RET_POTFIT_THRMO_210 if(nconfigs==210) return ThermoDynamics(POTFIT_THRMO_210); RET_POTFIT_THRMO_211
#define RET_POTFIT_THRMO_211 if(nconfigs==211) return ThermoDynamics(POTFIT_THRMO_211); RET_POTFIT_THRMO_212
#define RET_POTFIT_THRMO_212 if(nconfigs==212) return ThermoDynamics(POTFIT_THRMO_212); RET_POTFIT_THRMO_213
#define RET_POTFIT_THRMO_213 if(nconfigs==213) return ThermoDynamics(POTFIT_THRMO_213); RET_POTFIT_THRMO_214
#define RET_POTFIT_THRMO_214 if(nconfigs==214) return ThermoDynamics(POTFIT_THRMO_214); RET_POTFIT_THRMO_215
#define RET_POTFIT_THRMO_215 if(nconfigs==215) return ThermoDynamics(POTFIT_THRMO_215); RET_POTFIT_THRMO_216
#define RET_POTFIT_THRMO_216 if(nconfigs==216) return ThermoDynamics(POTFIT_THRMO_216); RET_POTFIT_THRMO_217
#define RET_POTFIT_THRMO_217 if(nconfigs==217) return ThermoDynamics(POTFIT_THRMO_217); RET_POTFIT_THRMO_218
#define RET_POTFIT_THRMO_218 if(nconfigs==218) return ThermoDynamics(POTFIT_THRMO_218); RET_POTFIT_THRMO_219
#define RET_POTFIT_THRMO_219 if(nconfigs==219) return ThermoDynamics(POTFIT_THRMO_219); RET_POTFIT_THRMO_220
#define RET_POTFIT_THRMO_220 if(nconfigs==220) return ThermoDynamics(POTFIT_THRMO_220); RET_POTFIT_THRMO_221
#define RET_POTFIT_THRMO_221 if(nconfigs==221) return ThermoDynamics(POTFIT_THRMO_221); RET_POTFIT_THRMO_222
#define RET_POTFIT_THRMO_222 if(nconfigs==222) return ThermoDynamics(POTFIT_THRMO_222); RET_POTFIT_THRMO_223
#define RET_POTFIT_THRMO_223 if(nconfigs==223) return ThermoDynamics(POTFIT_THRMO_223); RET_POTFIT_THRMO_224
#define RET_POTFIT_THRMO_224 if(nconfigs==224) return ThermoDynamics(POTFIT_THRMO_224); RET_POTFIT_THRMO_225
#define RET_POTFIT_THRMO_225 if(nconfigs==225) return ThermoDynamics(POTFIT_THRMO_225); RET_POTFIT_THRMO_226
#define RET_POTFIT_THRMO_226 if(nconfigs==226) return ThermoDynamics(POTFIT_THRMO_226); RET_POTFIT_THRMO_227
#define RET_POTFIT_THRMO_227 if(nconfigs==227) return ThermoDynamics(POTFIT_THRMO_227); RET_POTFIT_THRMO_228
#define RET_POTFIT_THRMO_228 if(nconfigs==228) return ThermoDynamics(POTFIT_THRMO_228); RET_POTFIT_THRMO_229
#define RET_POTFIT_THRMO_229 if(nconfigs==229) return ThermoDynamics(POTFIT_THRMO_229); RET_POTFIT_THRMO_230
#define RET_POTFIT_THRMO_230 if(nconfigs==230) return ThermoDynamics(POTFIT_THRMO_230); RET_POTFIT_THRMO_231
#define RET_POTFIT_THRMO_231 if(nconfigs==231) return ThermoDynamics(POTFIT_THRMO_231); RET_POTFIT_THRMO_232
#define RET_POTFIT_THRMO_232 if(nconfigs==232) return ThermoDynamics(POTFIT_THRMO_232); RET_POTFIT_THRMO_233
#define RET_POTFIT_THRMO_233 if(nconfigs==233) return ThermoDynamics(POTFIT_THRMO_233); RET_POTFIT_THRMO_234
#define RET_POTFIT_THRMO_234 if(nconfigs==234) return ThermoDynamics(POTFIT_THRMO_234); RET_POTFIT_THRMO_235
#define RET_POTFIT_THRMO_235 if(nconfigs==235) return ThermoDynamics(POTFIT_THRMO_235); RET_POTFIT_THRMO_236
#define RET_POTFIT_THRMO_236 if(nconfigs==236) return ThermoDynamics(POTFIT_THRMO_236); RET_POTFIT_THRMO_237
#define RET_POTFIT_THRMO_237 if(nconfigs==237) return ThermoDynamics(POTFIT_THRMO_237); RET_POTFIT_THRMO_238
#define RET_POTFIT_THRMO_238 if(nconfigs==238) return ThermoDynamics(POTFIT_THRMO_238); RET_POTFIT_THRMO_239
#define RET_POTFIT_THRMO_239 if(nconfigs==239) return ThermoDynamics(POTFIT_THRMO_239); RET_POTFIT_THRMO_240
#define RET_POTFIT_THRMO_240 if(nconfigs==240) return ThermoDynamics(POTFIT_THRMO_240); RET_POTFIT_THRMO_241
#define RET_POTFIT_THRMO_241 if(nconfigs==241) return ThermoDynamics(POTFIT_THRMO_241); RET_POTFIT_THRMO_242
#define RET_POTFIT_THRMO_242 if(nconfigs==242) return ThermoDynamics(POTFIT_THRMO_242); RET_POTFIT_THRMO_243
#define RET_POTFIT_THRMO_243 if(nconfigs==243) return ThermoDynamics(POTFIT_THRMO_243); RET_POTFIT_THRMO_244
#define RET_POTFIT_THRMO_244 if(nconfigs==244) return ThermoDynamics(POTFIT_THRMO_244); RET_POTFIT_THRMO_245
#define RET_POTFIT_THRMO_245 if(nconfigs==245) return ThermoDynamics(POTFIT_THRMO_245); RET_POTFIT_THRMO_246
#define RET_POTFIT_THRMO_246 if(nconfigs==246) return ThermoDynamics(POTFIT_THRMO_246); RET_POTFIT_THRMO_247
#define RET_POTFIT_THRMO_247 if(nconfigs==247) return ThermoDynamics(POTFIT_THRMO_247); RET_POTFIT_THRMO_248
#define RET_POTFIT_THRMO_248 if(nconfigs==248) return ThermoDynamics(POTFIT_THRMO_248); RET_POTFIT_THRMO_249
#define RET_POTFIT_THRMO_249 if(nconfigs==249) return ThermoDynamics(POTFIT_THRMO_249); RET_POTFIT_THRMO_250
#define RET_POTFIT_THRMO_250 if(nconfigs==250) return ThermoDynamics(POTFIT_THRMO_250); RET_POTFIT_THRMO_251
#define RET_POTFIT_THRMO_251 if(nconfigs==251) return ThermoDynamics(POTFIT_THRMO_251); RET_POTFIT_THRMO_252
#define RET_POTFIT_THRMO_252 if(nconfigs==252) return ThermoDynamics(POTFIT_THRMO_252); RET_POTFIT_THRMO_253
#define RET_POTFIT_THRMO_253 if(nconfigs==253) return ThermoDynamics(POTFIT_THRMO_253); RET_POTFIT_THRMO_254
#define RET_POTFIT_THRMO_254 if(nconfigs==254) return ThermoDynamics(POTFIT_THRMO_254); RET_POTFIT_THRMO_255
#define RET_POTFIT_THRMO_255 if(nconfigs==255) return ThermoDynamics(POTFIT_THRMO_255); RET_POTFIT_THRMO_256
#define RET_POTFIT_THRMO_256 return ThermoDynamics(POTFIT_THRMO_256);
#define RET_POTFIT_THRMO RET_POTFIT_THRMO_0
