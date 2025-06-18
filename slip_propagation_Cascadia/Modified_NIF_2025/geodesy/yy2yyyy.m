function YYYY=yy2yyyy(YY,pivot)
%YY2YYYY   Converts two digit years to four digits.
%     YYYY=yy2yyyy(YY,PIVOT) converts the matrix of
%     two digit years YY two a matrix of four digit
%     years YYYY, using the two digit PIVOT as the
%     pivot year.  PIVOT defaults to 80 if omitted.

%-------------------------------------------------------------------------------
%   Record of revisions:
%
%   Date          Programmer            Description of Change
%   ====          ==========            =====================
%
%   Aug 20, 2001  Peter Cervelli        Standardized code
%   Unknown       Peter Cervelli		Original Code
%
%-------------------------------------------------------------------------------

%Resort to default pivot if necessary

    if nargin<2
        pivot=80;
    end

%Convert

    I=YY<pivot;
    J=YY>=pivot & YY < 100;

    YYYY=YY;
    YYYY(I)=YYYY(I)+2000;
    YYYY(J)=YYYY(J)+1900;