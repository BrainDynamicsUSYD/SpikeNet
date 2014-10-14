function lab=f_lab(vect,num_values,origin,cut)  

if origin==1
    vect2=[0 vect];
else
    vect2=vect;
end
range=max(vect2)-min(vect2);
interval=range/num_values;
min_unit=fix(interval*10^(ceil(-log10(interval))))/10^(ceil(-log10(interval)));

xsep=round((range+min_unit)/num_values/min_unit)*min_unit;
lab=round(min(vect2)/min_unit)*min_unit+(xsep:xsep:(range+min_unit));

if cut==1
    lab=lab(lab>=min(vect2) & lab<=max(vect2));
end