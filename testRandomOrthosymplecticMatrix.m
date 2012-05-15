function testRandomOrthosymplecticMatrix

reset(RandStream.getDefaultStream);

assertExceptionThrown(@() randomOrthosymplecticMatrix(3),'cbrpack:oddSize');

n=5;

for k=1:25
S=randomOrthosymplecticMatrix(2*n);
assertVectorsAlmostEqual(S'*S,eye(2*n));
assertVectorsAlmostEqual(S'*jay(2*n)*S,jay(2*n));
end

for k=1:25
S=randomLagrangianSubspace(2*n);
assertVectorsAlmostEqual(S'*jay(2*n)*S,zeros(n));
end
