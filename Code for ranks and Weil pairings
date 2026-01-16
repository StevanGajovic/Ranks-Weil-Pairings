// Elliptic curves y^2=x^3+a^2
for a in [2, 3, 15, 427, 13*19*23*43] do
E := EllipticCurve([0, a^2]);
print Rank(E);
end for;

// Elliptic curves y^2=x^3+a^3
for a in [1, 2, 37, 506] do
E := EllipticCurve([0, a^3]);
print Rank(E);
end for;

// Elliptic curves y^2=(x+k^2-k)(x+k^2-1)(x+k^2+k)
for k in [2, 3, 11, 43, 329] do
Q<x> := PolynomialRing(Rationals());
E := EllipticCurve((x+k^2-k)*(x+k^2-1)*(x+k^2+k));
print Rank(E);
end for;

//Hyperelliptic Curves used in the construction for genus 3-6 curves - checking their ranks

function GiveRankg2(L) //we use this function because RankBounds(J) did not work for the curve y^2=(x+49)(x+4)(x-1)(x-43)(x-51)
Q<x> := PolynomialRing(Rationals());
f := &* [(x-l) : l in L];
C := HyperellipticCurve(f);
J := Jacobian(C);
return MordellWeilGroup(J);
end function;

GiveRankg2([1,2,3,4,10]);
GiveRankg2([1,2,3,4,37]);
GiveRankg2([-49,-4,1,43,51]);
GiveRankg2([-49,1,2,28,51]);


function GiveRankg3(L)
Q<x> := PolynomialRing(Rationals());
f := &* [(x-l) : l in L];
C := HyperellipticCurve(f);
J := Jacobian(C);
return RankBounds(J);
end function;

GiveRankg3([1,2,3,4,5,6,7]);
GiveRankg3([1,2,3,4,6,7,11]);
GiveRankg3([1,2,3,4,5,14,21]);

// Checking that the Weil Pairing on the basis of J(C) is indeed -1 for y^2=x(x-a_1)…(x-a_{2g+1}), where L=[a_1,…,a_{2g+1}] and a_1..a_{2g+1} is not zero - sufficient to check when contructing examples of genus 4 and 6 from the curves given in the paper
function CheckWP(L)
Q<x> := PolynomialRing(Integers());
A := &* [k: k in L];
f2 := &* [(x+Integers()!(A/k)) : k in L];
C2 := HyperellipticCurve(f2);
d := Integers()!(Discriminant(f2));
p := 1;
while not (IsPrime(p) and p mod 4 eq 1 and d mod p ne 0) do
p := p + 1;
end while;
Fp := GF(p);          
Fp<x> := PolynomialRing(Fp);
f2p := &+[Coefficients(f2)[i] mod p * x^(i-1) : i in [1..#Coefficients(f2)]];
C2p := HyperellipticCurve(f2p);
J2 := Jacobian(C2p);
LG := [J2![x+Integers()!(A/k),0] : k in L];
for i in [1..#L-2] do
for j in [i+1..#L-1] do
if Fp!(WeilPairing(LG[i],LG[j],2)+1) ne 0 then
return 0;
end if;
end for;
end for;
return 1;
end function;

CheckWP([1,2,3,4,10]);
CheckWP([1,2,3,4,37]);
CheckWP([-49,-4,1,43,51]);
CheckWP([-49,1,2,28,51]);

CheckWP([1,2,3,4,5,6,7]);
CheckWP([1,2,3,4,6,7,11]);
CheckWP([1,2,3,4,5,14,21]);

//Checking that the Weil Pairing on the basis of J(C) is indeed -1 for y^2=(x-(a_1-a_{2g+1}))…(x-(a_{2g}-a_{2g+1})), where L=[a_1,…,a_{2g+1}] and a_1..a_{2g+1} is not zero, by moving a_1-a_{2g+1} to zero and using the model y^2=x(x-(a_2-a_1))…(x-(a_{2g}-a_1)), and then applying CheckWP([a_2-a_1,…,a_{2g}-a_2]) - sufficient to check when contructing examples of genus 5 curves of genus 3 given in the paper
function CheckWP2(L) 
Q<x> := PolynomialRing(Rationals());
K := [(k-L[1]): k in L];
K2 := Exclude(K,K[#L]);
K3 := Exclude(K2,K2[1]);
return CheckWP(K3);
end function;

CheckWP2([1,2,3,4,5,6,7]);
CheckWP2([1,2,3,4,6,7,11]);
CheckWP2([1,2,3,4,5,14,21]);
