---
permalink: /docs/quick-start-guide/
title: "Quickstart"
excerpt: "A quick start guide."
sidebar:
    nav: docs
---

## This is a test page

```matlab
totalEdge = uint32(sort([elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])],2));
[edge,i2,j] = myunique(totalEdge);
NT = size(elem,1);
elem2edge = uint32(reshape(j,NT,3));
i1(j(3*NT:-1:1)) = 3*NT:-1:1; 
i1 = i1';
k1 = ceil(i1/NT); 
k2 = ceil(i2/NT); 
t1 = i1 - NT*(k1-1);
t2 = i2 - NT*(k2-1);
ix = (i1 ~= i2); 
neighbor = uint32(accumarray([[t1(ix),k1(ix)];[t2,k2]],[t2(ix);t1],[NT 3]));
edge2elem = uint32([t1,t2,k1,k2]);
bdElem = t1(t1 == t2);
bdk1 = k1(t1 == t2);
bdEdge = [elem(bdElem(bdk1==1),[2 3]); elem(bdElem(bdk1==2),[3 1]);...
          elem(bdElem(bdk1==3),[1 2])];
bdEdge2elem = [bdElem(bdk1==1);bdElem(bdk1==2);bdElem(bdk1==3)];
T = struct('neighbor',neighbor,'elem2edge',elem2edge,'edge',edge,...
           'edge2elem',edge2elem,'bdElem',bdElem,'bdEdge',bdEdge,...
           'bdEdge2elem', bdEdge2elem);
```

{% highlight matlab linenos %}

totalEdge = uint32(sort([elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])],2));
[edge,i2,j] = myunique(totalEdge);
NT = size(elem,1);
elem2edge = uint32(reshape(j,NT,3));
i1(j(3*NT:-1:1)) = 3*NT:-1:1; 
i1 = i1';
k1 = ceil(i1/NT); 
k2 = ceil(i2/NT); 
t1 = i1 - NT*(k1-1);
t2 = i2 - NT*(k2-1);
ix = (i1 ~= i2); 
neighbor = uint32(accumarray([[t1(ix),k1(ix)];[t2,k2]],[t2(ix);t1],[NT 3]));
edge2elem = uint32([t1,t2,k1,k2]);
bdElem = t1(t1 == t2);
bdk1 = k1(t1 == t2);
bdEdge = [elem(bdElem(bdk1==1),[2 3]); elem(bdElem(bdk1==2),[3 1]);...
          elem(bdElem(bdk1==3),[1 2])];
bdEdge2elem = [bdElem(bdk1==1);bdElem(bdk1==2);bdElem(bdk1==3)];
T = struct('neighbor',neighbor,'elem2edge',elem2edge,'edge',edge,...
           'edge2elem',edge2elem,'bdElem',bdElem,'bdEdge',bdEdge,...
           'bdEdge2elem', bdEdge2elem);

{% endhighlight %}