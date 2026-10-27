loadPackage "SymbolicPowers";
R = QQ[x,y];
I = ideal(x^3,x*y);
S = symbolicPower(I,2,UseMinimalPrimes=>true);
print("CHECK|minimal_square_is_x2|" | toString(S==ideal(x^2)));
print("CHECK|actual_minimal_containment|" | toString isSubset(S,I));
print("CHECK|reported_minimal_containment|" | toString isSymbPowerContainedinPower(I,2,1,UseMinimalPrimes=>true));
print("CHECK|reported_minimal_reverse_containment|" | toString containmentProblem(I,2,InSymbolic=>true,UseMinimalPrimes=>true));

R = QQ[x,y,z];
I = ideal(x*y,x*z,y*z);
print("CHECK|zeroth_power_is_unit|" | toString(symbolicPower(I,0)==ideal(1_R)));
print("CHECK|nonradical_join_returns_x3|" | toString(symbolicPowerJoin(ideal(x^2),2)==ideal(x^3)));

R = QQ[x,y,z,w];
J = ideal((x+y)*y,(x+y)*z,y*z);
I = intersect(J,(ideal vars R)^3);
print("CHECK|embedded_default_symbolic_square_is_ordinary|" | toString(symbolicPower(I,2)==I^2));
print("CHECK|old_associated_heights|" | toString assPrimesHeight I);
print("CHECK|new_square_associated_heights|" | toString assPrimesHeight(I^2));

B = QQ[a,b,c];
Q = B/ideal(a*b-c^2);
P = ideal(Q_0,Q_2);
print("CHECK|quadric_ambient_is_quotient|" | toString(ring P === Q));
print("CHECK|quadric_quotient_symbolic_square_is_a|" | toString(symbolicPower(P,2)==ideal(Q_0)));
print "CORRECTNESS_COMPLETE";
exit 0;
