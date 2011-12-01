function testRandomOrthosymplecticMatrix

reset(RandStream.getDefaultStream);

assertExceptionThrown(@() randomOrthosymplecticMatrix(3),'cbrpack:oddSize');

n=5;

for k=1:100
S=randomOrthosymplecticMatrix(2*n);
assertVectorsAlmostEqual(S*S',eye(2*n));
assertVectorsAlmostEqual(S*matgic.jay(2*n)*S',matgic.jay(2*n));
end
