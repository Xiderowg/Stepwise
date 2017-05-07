clc
clear all
y=importdata('预报对象.txt');%预报对象
x=importdata('影响因子.txt');%预报因子

M=[x y];%合并预报因子和预报对象
n=size(M,1);%数据总数
XNUM=size(x,2);%预报因子数目
TOTALNUM=size(M,2);%因子和对象总数
MM=roundn(mean(M),-2);%平均值矩阵
R=corrcoef(M);%相关系数矩阵
COVM=cov(M);%协方差
flag=0;%循环次数
INDEXG=[];%引入因子数组
while 1
    flag=flag+1;
    V=[];%所有因子方差贡献量数组
    for i=1:XNUM
        V=[V power(R(i,TOTALNUM),2)/R(i,i)];
    end
    
    if flag>3 %第三次引入因子后需要进行剔除检验
        %[VM,INDEX]=min(V);
        VM=V(1);
        for tt=1:XNUM
            if V(tt)<=VM && isIn(INDEXG,tt)
                INDEX=tt;
                VM=V(tt);
            end
        end
        str=['Check:x' num2str(INDEX)];
        disp(str);
        %VM为以选中因子中最小的方差贡献值，INDEX为贡献值的下标
        F2=VM*(n-length(INDEXG)-1)/R(TOTALNUM,TOTALNUM);
        if F2<=5.0
            K=INDEX;
            %删除因子
            for tt=1:size(INDEXG,2)
                if INDEXG(tt)==K
                    INDEXG(tt)=[];
                    str=['Delete:x' num2str(tt)];
                    disp(str);
                end
            end
        else
            %[VM,INDEX]=max(V);
            VM=0;
            for tt=1:XNUM
                if V(tt)>=VM && isIn(INDEXG,tt)==false
                    INDEX=tt;
                    VM=V(tt);
                end
            end
         %VM为以未选中因子中最大的方差贡献值，INDEX为贡献值的下标
            F1=VM*(n-length(INDEXG)-2)/(R(TOTALNUM,TOTALNUM)-VM);
            if F1<=5.0 %既不能剔除也不能加入
                break;
            else
                K=INDEX;
                INDEXG=[INDEXG K];
                str=['Introduce:x' num2str(K)];
                disp(str);
            end
        end
    else
        %[VM,INDEX]=max(V);
        VM=0;
        for tt=1:XNUM
          if (V(tt)>=VM) && (isIn(INDEXG,tt)==false)
              INDEX=tt;
              VM=V(tt);
          end
        end
        %VM为以未选中因子中最大的方差贡献值，INDEX为贡献值的下标
        F1=VM*(n-length(INDEXG)-2)/(R(TOTALNUM,TOTALNUM)-VM);
        if F1<=5.0 %既不能剔除也不能VN入
            break;
        else
            K=INDEX;
            INDEXG=[INDEXG K];
            str=['Introduce:x' num2str(K)];
            disp(str);
        end
    end
    RTMP=R;%上一次变化的矩阵
    for x=1:TOTALNUM
        for y=1:TOTALNUM
            if (x~=K) && (y~=K)
                R(x,y)=RTMP(x,y)-RTMP(x,K)*RTMP(K,y)/RTMP(K,K);
            else if x==y
                    R(K,K)=1/RTMP(K,K);
                else if y==K
                        R(x,K)=-RTMP(x,K)/RTMP(K,K);
                    else
                        R(K,y)=RTMP(K,y)/RTMP(K,K);
                    end
                end
            end
        end
    end
    
end
disp('End of regression');
streqa=['y='];
b=[MM(size(MM,2))];
for ii=1:size(INDEXG,2)
   str=['b_',num2str(ii),'='];
   b=[b 0];
   disp(str);
   disp(R(INDEXG(ii),TOTALNUM));
   b(1)=b(1)-R(INDEXG(ii),TOTALNUM)*sqrt(COVM(TOTALNUM,TOTALNUM))/sqrt(COVM(INDEXG(ii),INDEXG(ii)))*MM(INDEXG(ii));
   b(ii+1)=R(INDEXG(ii),TOTALNUM)*sqrt(COVM(TOTALNUM,TOTALNUM))/sqrt(COVM(INDEXG(ii),INDEXG(ii)));
   streqa=[streqa,'(',num2str(b(ii+1)),')','*x',num2str(INDEXG(ii)),'+'];
end
streqa=[streqa,'(',num2str(b(1)),')'];
disp('Final equation:')
disp(streqa);