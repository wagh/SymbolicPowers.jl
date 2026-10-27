-- Appended to the common benchmark definitions; constructors return fresh ideals.
RF5 = (ZZ/5)[k_0..k_2, MonomialOrder=>GRevLex];
RF101 = (ZZ/101)[l_0..l_2, MonomialOrder=>GRevLex];
pointCoordinates = {{1,0,0},{0,1,0},{0,0,1},{1,1,1},{1,2,3},{3,1,-2},{2,3,1},{1,-1,2}};
pointIdeal = (R,v) -> (
    j := first select(0..2, i -> v#i != 0);
    ideal apply(select(0..2, i -> i != j), i -> v#j*R_i-v#i*R_j)
    );
eightPoints = () -> ideal gens intersect apply(pointCoordinates, v -> pointIdeal(R3,v));
fatPoints = () -> ideal gens intersect apply(0..3, i -> (pointIdeal(R3,pointCoordinates#i))^(if i==0 then 2 else 1));
mixedIdeal = () -> ideal gens intersect(ideal(R4_0^2+R4_1^2-R4_2),ideal(R4_0-1,R4_1-R4_2^2,R4_3));
extraCases = {
    {"points8_sp2",eightPoints,I -> symbolicPower(I,2),"ideal"},
    {"points8_sp4",eightPoints,I -> symbolicPower(I,4),"ideal"},
    {"fatpoints4_sp2",fatPoints,I -> symbolicPower(I,2),"ideal"},
    {"points8_reg3",eightPoints,I -> asymptoticRegularity(I,SampleSize=>3),"value"},
    {"fermat_ci3",() -> fermat R3,I -> symbolicPower(I,3,CIPrimes=>true),"ideal"},
    {"mixed_sp2",mixedIdeal,I -> symbolicPower(I,2),"ideal"},
    {"frobenius_linear101_sp2",() -> ideal(RF101_0+RF101_1,RF101_2),I -> symbPowerPrimePosChar(I,2),"ideal"},
    {"frobenius_point5_sp2",() -> ideal(RF5_0-1,RF5_1,RF5_2+RF5_0),I -> symbPowerPrimePosChar(I,2),"ideal"}
    };
monoIds = {"triangle_sp6","triangle_sp12","embedded_sp3","ci_monomial_sp8",
    "contain_triangle8","equal_triangle4","wald_triangle","wald_cycle5","defect_triangle2",
    "resurgence_triangle5","reg_triangle4","konig_complete6","konig_hyper6",
    "packed_cycle4","packed_cycle5","packed_triangle_ambient8","squarefree_m5","squarefree_m6"};
monoCases = apply(select(cases, spec -> member(spec#0,monoIds)), spec -> {
    "mono_"|spec#0, () -> monomialIdeal((spec#1)()), spec#2, spec#3});
cases = join(select(cases,spec -> not member(spec#0,{"contain_triangle24","squarefree_m7","squarefree_m10"})),extraCases,monoCases);
