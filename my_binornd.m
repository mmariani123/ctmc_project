%Draws random number from a binomial distribution, using a normal
%approximation for appropriate values using the condition
%"abs(1./sqrt(n)*(sqrt((1-p)./p)-sqrt(p./(1-p))))<0.3"
function R=my_binornd(varargin)
%function R=my_binornd(varargin)
assert(nargin>=2)%at least two input parameters are necessary
n=varargin{1};
p=varargin{2};
if(length(varargin)>2)
  more_args=cell2mat(varargin(3:end));
else
  %more_args=double.empty();
  %more_args=double(0);
  more_args=[];
end

if(numel(n)==1)
  n=n*ones(size(p));
elseif(numel(p)==1),p=p*ones(size(n));end

% criteria=abs(1./sqrt(n).*(sqrt((1-p)./p)-sqrt(p./(1-p))))<0.3 & n>5;%criteria for normal approximation (Box, Hunter and Hunter (1978). Statistics for experimenters. Wiley. p. 130)
criteria=n.*p>5 & n.*(1-p)>5;%criteria for normal approximation
n1=n(criteria);
n2=n(~criteria&n>0);
p1=p(criteria);
p2=p(~criteria&n>0);

if(isempty(more_args))
    if(isempty(n1)),R1=[];
    else R1=normrnd(n1.*p1,sqrt(n1.*p1.*(1-p1)));
        while(nnz(R1<0)>0)%re-draw negative numbers
            a=R1<0;
            R1(a)=normrnd(n1(a).*p1(a),sqrt(n1(a).*p1(a).*(1-p1(a))));
        end
    end%normal approximation

    if(isempty(n2)),R2=[];
    else R2=binornd(n2,p2);
    end%binomial draw
else
    if(isempty(n1)),R1=[];
    else R1=normrnd(n1.*p1,sqrt(n1.*p1.*(1-p1)),more_args);
        while(nnz(R1<0)>0)%re-draw negative numbers
            a=R1<0;
            R1(a)=normrnd(n1(a).*p1(a),sqrt(n1(a).*p1(a).*(1-p1(a))),more_args);
        end
    end%normal approximation

    if(isempty(n2)),R2=[];
    else R2=binornd(n2,p2,more_args);
    end%binomial draw
end

R1=round(R1);
R=zeros(size(n));
if(length(criteria)==1)
    if(criteria),R=R1;
    else R=R2;
    end
else
    R(criteria)=R1;
    R(~criteria&n>0)=R2;
    %note - this way the remaining values, where n=0, are set to 0
end

%to see the approximation, run this line:
% figure;hold on;n=1:100;for p=[0.1 0.2 0.3],plot(n,abs(1./sqrt(n).*(sqrt((1-p)./p)-sqrt(p./(1-p)))),'Color',rand(1,3));end,legend({'p=0.1','p=0.2','p=0.3'}),plot([1 length(n)],[0.3 0.3],'--')


%my_binornd against binornd
%mb=my_binornd(repmat([2 100],1e6,1),repmat([0.5 0.1],1e6,1));b=binornd(repmat([2 100],1e6,1),repmat([0.5 0.1],1e6,1));
%figure;subplot(2,2,1),hist(mb(:,1)',-1:5),;subplot(2,2,2),hist(b(:,1)',-1:5),subplot(2,2,3),hist(mb(:,2)',-1:40),;subplot(2,2,4),hist(b(:,2)',-1:40)
