* FLAT in the cloud
* Land allocation tool prepared for HubZero
* Coded by Jingyu Song
* 02/12/2016
* Updated 05/21/2016 to include projected fraction calculations
* This method is based on Papke and Wooldridge (1993)
* We adjusted the Quasi-maximum likelihood estimator they proposed to our problem

$offdigit
option limrow=30,limcol=0 ;
OPTION Decimals=8;

********************************************************************************
** Load data. Data files are generated based on the dataset provided to FLAT
********************************************************************************
SET state States
/

$offlisting
$include states.csv
$onlisting

/;
SET pixel  pixels
/
$offlisting
$include pixels.csv
$onlisting
/;

SET k variable names
/

$offlisting
$include names.csv
$onlisting

/;
DISPLAY k;

SET crop crop names
/

$offlisting
$include cropnames.csv
$onlisting

/;
DISPLAY crop;

PARAMETER statearea(state) state level total area information
/
$offlisting
$ondelim
$include statelevelareainfo.csv
$offdelim
$onlisting
/;

TABLE data(pixel,state,k) independent variables
$offlisting
$ondelim
$include data.csv
$offdelim
$onlisting

TABLE areainfo(state,crop) state level crop area info
$offlisting
$ondelim
$include statelevelcroparea.csv
$offdelim
$onlisting

set header /lon,lat,pixelarea/;

TABLE pixelarea(pixel,state,header) pixel area calculated based on lon and lat
$offlisting
$ondelim
$include pixelarea.csv
$offdelim
$onlisting

********************************************************************************
** Model calculations
********************************************************************************
PARAMETER cropfraction(state,crop);
cropfraction(state,crop) = areainfo(state,crop)/statearea(state);


PARAMETER totalcropfraction(state);
totalcropfraction(state) = sum(crop,cropfraction(state,crop));
DISPLAY totalcropfraction;

PARAMETER totalcroparea(state);
totalcroparea(state) = sum(pixel,pixelarea(pixel,state,'pixelarea'));


SET map(pixel,state) map pixels with states;
map(pixel,state) = yes $ sum(k,data(pixel,state,k));

VARIABLES  theta(crop,k)  coefficients on each land property variable
           z              likelihood function to be maximized ;

EQUATION obj objective function of likelihood to maximize ;

ALIAS(crop,crop1) ;

********************************************************************************
** Quasi-maximum likelihood set up
********************************************************************************

obj..z =e=
           sum(state$(totalcroparea(state)),
           (1-totalcropfraction(state))
           * log(1-sum(pixel$map(pixel,state),sum(crop, exp(sum(k,theta(crop,k)*data(pixel,state,k)))
           /(1+sum(crop1,exp(sum(k,theta(crop1,k)*data(pixel,state,k)))))*pixelarea(pixel,state,'pixelarea')))/statearea(state))

           + sum(crop,
           cropfraction(state,crop)
           * log(sum(pixel$map(pixel,state),exp(sum(k,theta(crop,k)*data(pixel,state,k)))
           /(1+sum(crop1,exp(sum(k,theta(crop1,k)*data(pixel,state,k)))))*pixelarea(pixel,state,'pixelarea'))/statearea(state))
           ))
           ;
MODEL likelihood/obj/;

SOLVE likelihood using nlp maximizing z;
DISPLAY z.l,theta.l;



********************************************************************************
** Predicted fractions/results output
********************************************************************************
PARAMETER fractions(pixel,state,crop) fractions at pixel level;
fractions(pixel,state,crop)$(totalcroparea(state) and map(pixel,state)) = exp(sum(k,theta.l(crop,k)
           *data(pixel,state,k)$map(pixel,state)))
           /(1+sum(crop1,exp(sum(k,theta.l(crop1,k)*data(pixel,state,k)))));
*DISPLAY fractions;

PARAMETER longitude(pixel);
longitude(pixel) = sum(state,pixelarea(pixel,state,'lon'));
*DISPLAY longitude;

PARAMETER latitude(pixel);
latitude(pixel) = sum(state,pixelarea(pixel,state,'lat'));
*DISPLAY latitude;

********************************************************************************
** Scale predicted fractions
********************************************************************************

* By scaling the predicted fractions,the predicted state total matches FAO state total
ALIAS (pixel,pixel2);
PARAMETER scaledfractions(pixel,state,crop) scaled fractions at pixel level;
scaledfractions(pixel,state,crop)$fractions(pixel,state,crop) = fractions(pixel,state,crop)*
           areainfo(state,crop)/sum(pixel2$map(pixel2,state),fractions(pixel2,state,crop)*pixelarea(pixel2,state,'pixelarea'));

*DISPLAY scaledfractions;


********************************************************************************
** Predicted fractions are exported to file finalresults.dat
********************************************************************************
file finalresults/finalresults.dat/ ;
put finalresults ;



put ' Pixel' @15;
put 'State' @33;
put 'lon' @46;
put 'lat' @54;

loop(crop$(ord(crop) le card(crop)),
     put 'Crop','    ', 'Fraction','     ' ;
);
put /;

loop ((pixel,state),
if (map(pixel,state),

     put pixel.tl:11,' ',state.tl:15,' ';
     put longitude(pixel):11:5;
     put latitude(pixel):11:5;


loop(crop,

if ((ord(crop) eq 1),

     put '  ',crop.tl:7,' ';
     put scaledfractions(pixel,state,crop):11:8;

else
     put '  ',crop.tl:7,' ';
     put scaledfractions(pixel,state,crop):11:8;

);

);
put /;

);
 );


********************************************************************************
** Calculate the hessian
********************************************************************************
VARIABLE zhessian;
EQUATION objhessian objective function of likelihood to maximize
         nhessian(state)    hessian for each of the selected state;
objhessian..zhessian =e= 0   ;
nhessian(state) $ totalcroparea(state)..
           (1-totalcropfraction(state))
           * log(1-sum(pixel$map(pixel,state),sum(crop, exp(sum(k,theta(crop,k)*data(pixel,state,k)))
           /(1+sum(crop1,exp(sum(k,theta(crop1,k)*data(pixel,state,k)))))*pixelarea(pixel,state,'pixelarea')))/statearea(state))

           + sum(crop,
           cropfraction(state,crop)
           * log(sum(pixel$map(pixel,state),exp(sum(k,theta(crop,k)*data(pixel,state,k)))
           /(1+sum(crop1,exp(sum(k,theta(crop1,k)*data(pixel,state,k)))))*pixelarea(pixel,state,'pixelarea'))/totalcroparea(state))
           )
           =n= 0;


MODEL likelihoodhessian/objhessian,nhessian/;
likelihoodhessian.optfile = 1 ;
option nlp=conopt ;
solve likelihoodhessian using nlp maximizing zhessian ;



set dummy / e1*e10000,x1*x10000 / ;

SET i(*),j(*) ;
PARAMETER a(*,*);
Variable x(*);
PARAMETER h(*,*,*) ;

Execute_Load 'hessian',i,j,a,x.l,h ;
*display i,j,a,x.l,h ;

********************************************************************************
** Calculate and output the covariance matrix, estimates and se
********************************************************************************

SET index variables
/

$offlisting
$include variables.csv
$onlisting

/;


PARAMETER hessian(*,index,index);
ALIAS (index,indexi,indexii,indexiii);
LOOP(indexii,
hessian(i,index,indexii)=h(i,index,indexii);
hessian(i,indexii,index)=hessian(i,index,indexii);
);
*DISPLAY hessian;



PARAMETER hessian_average(index,index);
hessian_average(index,indexii) = 1/(card(i)-1)*sum(i,hessian(i,index,indexii));
*DISPLAY hessian_average;

PARAMETER hessian_inverse(index,index);
execute_unload 'inverse_hessian_average.gdx',index,hessian_average;
execute 'invert inverse_hessian_average.gdx index hessian_average inversehessianmaize.gdx hessian_inverse ';
execute_load 'inversehessianmaize.gdx',hessian_inverse;
*DISPLAY hessian_inverse;


PARAMETER B(index,indexii) outer product of the gradient;
a('e1',index)=0;
B(index,indexii) = 1/(card(i)-1)*sum(i,a(i,index)*a(i,indexii));
*DISPLAY B;

PARAMETER Covariance(index,index);
Covariance(index,indexi) = sum(indexiii,hessian_inverse(index,indexiii)*sum(indexii,B(indexiii,indexii)*
                           hessian_inverse(indexii,indexi)))/(card(i)-1);
*DISPLAY Covariance;

** append 'real' variable names to estimation results

SET names variable names
/

$offlisting
$include names.csv
$onlisting

/;
*DISPLAY names;

SET cropnames crop names
/

$offlisting
$include cropnames.csv
$onlisting

/;
*DISPLAY cropnames;

********************************************************************************
** the covariance matrix is exported to file covariances.dat
********************************************************************************
file covariances/covariances.dat/;
put covariances;
set q /1*10/;
loop(q$(card(index) gt (ord(q) - 1)*10),
* Column headers
  put '               ';   put '          ';
  loop(indexi$(ord(indexi) gt 10*(ord(q)-1) and ord(indexi) le 10*ord(q)),
      loop((cropnames,names)$((ord(cropnames)-1)*card(names)+ord(names) eq ord(indexi)),
        put cropnames.tl:14 );
    ) ;
    put / ;
  put '               ';   put '          ';
  loop(indexi$(ord(indexi) gt 10*(ord(q)-1) and ord(indexi) le 10*ord(q)),
      loop((cropnames,names)$((ord(cropnames)-1)*card(names)+ord(names) eq ord(indexi)),
        put names.tl:14 );
    ) ;
    put / ;

* Each row of data ;
loop(index,
  loop((cropnames,names)$((ord(cropnames)-1)*card(names)+ord(names) eq ord(index)),
    put cropnames.tl:10,names.tl:10 );
    loop(indexi$(ord(indexi) gt 10*(ord(q)-1) and ord(indexi) le 10*ord(q)),
      put Covariance(index,indexi):14:8 ;
    ) ;
    put / ;
  ) ;
put /;
) ;

PARAMETER se(index) standard error;
se(index) = sqrt(Covariance(index,index));
DISPLAY se;

PARAMETER t(index) t value;
t(index) = x.l(index)/se(index);

********************************************************************************
** the coefficient estimates and se are exported to file estimates.dat
********************************************************************************
file estimates/estimates.dat/;
put estimates;

put 'Crop' @13;
put 'Variable' @28;
put 'Coefficient' @45;
put 'Standard Error';

loop ((cropnames,names),
    put / cropnames.tl:11,'  ', names.tl:11;
);
loop ((index),
     put #(1+ord(index)),@25, x.l(index):14:8,'       '
     put ' ' se(index):12:8 ' ';
);
