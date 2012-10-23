%%tests
p=[1 2; 1 2; 1 2; 1 2];
m=[1 2 3 0;1 2 3 4]';
v=[2 3 4 3; 2 3 4 5]';
gmm=gmmCreate([1],[1 1],[.3 .7],m,v);
A=gmmLikelihood(p,gmm);

B=[mvnpdf(p(:,1)',[1 2 3 0],[2 3 4 3])*.3+mvnpdf(p(:,1)',[1 2 3 4],[2 3 4 5])*.7,... 
mvnpdf(p(:,2)',[1 2 3 0],[2 3 4 3])*.3+mvnpdf(p(:,2)',[1 2 3 4],[2 3 4 5])*.7];

assert(approxeq(A,B, 1e-15));


gmm=gmmCreate([.5;.5],[1 2],[1 1],m,v);
A=gmmLikelihood(p,gmm)*2;

B=[mvnpdf(p(:,1)',[1 2 3 0],[2 3 4 3]) mvnpdf(p(:,1)',[1 2 3 4],[2 3 4 5]) ,... 
mvnpdf(p(:,2)',[1 2 3 0],[2 3 4 3]) mvnpdf(p(:,2)',[1 2 3 4],[2 3 4 5]) ];

assert(approxeq(A,B, 1e-15));