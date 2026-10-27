-- SymbolicPowers 2.0 on Macaulay2: fresh ideals, warmed code, wall-clock time.
loadPackage "SymbolicPowers";
print("META|Macaulay2|" | toString version#"VERSION");
print("META|SymbolicPowers|" | toString SymbolicPowers.Options.Version);

R3 = QQ[x,y,z, MonomialOrder => GRevLex];
RW = QQ[u,v,w, Degrees => {3,4,5}, MonomialOrder => GRevLex];
R4 = QQ[a,b,c,d, MonomialOrder => GRevLex];
R5 = QQ[e_0..e_4, MonomialOrder => GRevLex];
R6 = QQ[f_0..f_5, MonomialOrder => GRevLex];
R7 = QQ[g_0..g_6, MonomialOrder => GRevLex];
R8 = QQ[h_0..h_7, MonomialOrder => GRevLex];
R10 = QQ[j_0..j_9, MonomialOrder => GRevLex];
RF = (ZZ/3)[p,q,r, MonomialOrder => GRevLex];

triangle = R -> ideal(R_0*R_1, R_0*R_2, R_1*R_2);
fermat = R -> (
    xx := R_0; yy := R_1; zz := R_2;
    ideal(xx*(yy^3-zz^3), yy*(zz^3-xx^3), zz*(xx^3-yy^3))
    );
curve = R -> (
    xx := R_0; yy := R_1; zz := R_2;
    ideal(yy^2-xx*zz, xx^2*yy-zz^2, xx^3-yy*zz)
    );
twisted = R -> ideal(R_0*R_2-R_1^2, R_0*R_3-R_1*R_2, R_1*R_3-R_2^2);
embedded = R -> ideal(R_0^3, R_0^2*R_1^2, R_0*R_2^3, R_1^4, R_1^2*R_2);
cycle = R -> ideal apply(0..numgens(R)-1, k -> R_k*R_((k+1)%numgens(R)));
hyper = R -> ideal apply(subsets(numgens R,3), s -> product apply(s,k -> R_k));

-- {identifier, ideal constructor, operation, result kind}
cases = {
    {"triangle_sp6", () -> triangle R3, I -> symbolicPower(I,6), "ideal"},
    {"triangle_sp12", () -> triangle R3, I -> symbolicPower(I,12), "ideal"},
    {"fermat_sp2", () -> fermat R3, I -> symbolicPower(I,2), "ideal"},
    {"fermat_sp3", () -> fermat R3, I -> symbolicPower(I,3), "ideal"},
    {"curve_sp3", () -> curve R3, I -> symbolicPower(I,3), "ideal"},
    {"curve_weighted_sp3", () -> curve RW, I -> symbolicPower(I,3), "ideal"},
    {"twisted_sp3", () -> twisted R4, I -> symbolicPower(I,3), "ideal"},
    {"embedded_sp3", () -> embedded R3, I -> symbolicPower(I,3), "ideal"},
    {"embedded_min_sp3", () -> embedded R3, I -> symbolicPower(I,3,UseMinimalPrimes=>true), "ideal"},
    {"ci_monomial_sp8", () -> ideal(R3_0^2,R3_1^2), I -> symbolicPower(I,8), "ideal"},
    {"affine_sp3", () -> ideal(R4_1*R4_2+1,R4_0*R4_2*R4_3), I -> symbolicPower(I,3), "ideal"},
    {"join_triangle4", () -> triangle R3, I -> symbolicPowerJoin(I,4), "ideal"},
    {"join_triangle4_core", () -> triangle R3, I -> symbolicPowerJoin(I,4), "ideal"},
    {"join_twisted3_core", () -> twisted R4, I -> symbolicPowerJoin(I,3), "ideal"},
    {"frobenius_curve3", () -> curve RF, I -> symbPowerPrimePosChar(I,3), "ideal"},
    {"frobenius_curve4", () -> curve RF, I -> symbPowerPrimePosChar(I,4), "ideal"},
    {"contain_triangle8", () -> triangle R3, I -> containmentProblem(I,8), "value"},
    {"contain_triangle24", () -> triangle R3, I -> containmentProblem(I,24), "value"},
    {"contain_fermat2", () -> fermat R3, I -> containmentProblem(I,2), "value"},
    {"equal_triangle4", () -> triangle R3, I -> isSymbolicEqualOrdinary(I,4), "value"},
    {"height_fermat3", () -> fermat R3, I -> bigHeight(I), "value"},
    {"wald_triangle", () -> triangle R3, I -> waldschmidt(I), "value"},
    {"wald_cycle5", () -> cycle R5, I -> waldschmidt(I), "value"},
    {"defect_triangle2", () -> triangle R3, I -> symbolicDefect(I,2), "value"},
    {"resurgence_triangle5", () -> triangle R3, I -> lowerBoundResurgence(I,SampleSize=>5,UseWaldschmidt=>true), "value"},
    {"reg_triangle4", () -> triangle R3, I -> asymptoticRegularity(I,SampleSize=>4), "value"},
    {"konig_complete6", () -> ideal vars R6, I -> isKonig(I), "value"},
    {"konig_hyper6", () -> hyper R6, I -> isKonig(I), "value"},
    {"packed_cycle4", () -> cycle R4, I -> isPacked(I), "value"},
    {"packed_cycle5", () -> cycle R5, I -> isPacked(I), "value"},
    {"packed_triangle_ambient8", () -> triangle R8, I -> isPacked(I), "value"},
    {"squarefree_m5", () -> ideal vars R5, I -> squarefreeInCodim(I), "monomials"},
    {"squarefree_m6", () -> ideal vars R6, I -> squarefreeInCodim(I), "monomials"},
    {"squarefree_m7", () -> ideal vars R7, I -> squarefreeInCodim(I), "monomials"},
    {"squarefree_m10", () -> ideal vars R10, I -> squarefreeInCodim(I), "monomials"}
    };

emitIdeal = (id, I) -> (
    print("IDEAL|" | id | "|" | toString numgens ring I | "|" | toString char ring I);
    scan(toList pairs flatten entries gens I, pair -> (
        idx := pair#0; ff := pair#1;
        print("POLY|" | id | "|" | toString idx);
        scan(terms ff, tt -> print("TERM|" | id | "|" | toString idx | "|" |
            toString leadCoefficient tt | "|" | toString first exponents tt));
        ));
    );

runCase = spec -> (
    id := spec#0; makeIdeal := spec#1; operation := spec#2; kind := spec#3;
    print("START|" | id);
    flush stdio;
    warmInput := makeIdeal(); warmAnswer := operation(warmInput);
    answer := null;
    samples := for trial from 1 to 5 list (
        inputIdeal := makeIdeal();
        measurement := elapsedTiming operation(inputIdeal);
        answer = measurement#1;
        measurement#0
        );
    ordered := sort samples;
    print("BENCH|" | id | "|" | toString ordered#2 | "|" | toString first ordered | "|" | toString last ordered | "|" | toString samples);
    if kind == "ideal" then emitIdeal(id,answer) else (
        if kind == "monomials" then emitIdeal(id,ideal answer) else print("VALUE|" | id | "|" | toString answer)
        );
    flush stdio;
    );

load "benchmark/comparison/extra-cases.m2";

if getenv "BENCH_DEFINITIONS_ONLY" =!= "1" then (
    selectedCase := getenv "BENCH_CASE";
    selectedCases := if selectedCase =!= "" and selectedCase =!= null then select(cases, spec -> spec#0 == selectedCase) else select(cases, spec -> not member(spec#0,{"contain_triangle24","squarefree_m7","squarefree_m10"}));
    failures := {};
    scan(selectedCases, spec -> try runCase(spec) else (print("ERROR|" | spec#0); failures=append(failures,spec#0)));
    if #failures > 0 then exit 1;
    print "BENCHMARKS_COMPLETE";
    exit 0;
    );
