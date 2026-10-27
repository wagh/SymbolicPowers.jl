loadPackage "SymbolicPowers";
R = QQ[x,y,Degrees=>{2,3}];
I = ideal(x^3,y^2);
print("WEIGHTED|same_ideal|" | toString(ideal monomialIdeal I == I));
print("WEIGHTED|ordinary_Ideal_waldschmidt|" | toString waldschmidt I);
print("WEIGHTED|MonomialIdeal_waldschmidt|" | toString waldschmidt monomialIdeal I);
print("WEIGHTED|expected_with_ring_weights|6");
exit 0;
