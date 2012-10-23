%unit tests for the pfile io functions

%%Test the reads
%old style
testF='small.pfile';
F=pfread(testF,0);
%new style
rf=pfopen(testF,'r');
Fnew=pfread(rf,0);
assert(all(F(:)==Fnew(:)))

[F,~,D]=pfread(testF,[0 3]);


%%Test the writes
[R,U,f,L,T] = pfinfo(testF,-1);
pfid=pfopen('writeTry.pfile','w');

bigF=[];
for i=0:U-1
    F=pfread(rf,i)';
    bigF=[bigF F];
    pfwrite(pfid,F);
end
pfclose(pfid);

rf2=pfopen('writeTry.pfile','r');
[R2,U2,f2,L2,T2] = pfinfo('writeTry.pfile',-1);
bigFAgain=[];
for i=0:U-1
    F=pfread(rf2,i)';
    bigFAgain=[bigFAgain F];
    if i ==98
        1;
    end
end
pfclose(rf2);
assert(all(bigFAgain(:)==bigF(:)));

%%test multiple reads and multiple writes
[F1,~,I1]=pfread(rf,1);
[F2,~,I2]=pfread(rf,2);
rf2=pfopen('writeTry2.pfile','w');
pfwrite(rf2,[F1' F2'],[I1' I2']);
pfclose(rf2);
rf2=pfopen('writeTry2.pfile','r');
[F3,~,I3]=pfread(rf2,[0 1]);
pfclose(rf2);
origF3=[F1;F2];
assert(all(origF3(:)==F3(:)));
origI3=[I1;I2];
origI3(:,1)=origI3(:,1)-1;
assert(all(origI3(:)==I3(:)));
pfclose(rf);
disp('tests passed')