%let path = U:\Consulting\KEL\Fox\Samuel;
libname sam "&path.\input";

data sheet1; set sam.sheet1; id = patient; run;
data sheet2; set sam.sheet2; run;
data sheet3; set sam.sheet3; run;





***** ANOVA Analysis;
%macro aov(ds,mod_ty,mod,research_q,iv,dv,delivery);

data ds; 
set &ds; 
%if "&ds" = "sheet3" %then %do;
where delivery = "&delivery."; 
%end;
run;


proc mixed data=ds;
		class &iv.;
		model &dv. = &iv.  / solution;
		lsmeans &iv.  / pdiff=all;
		ods output lsmeans=lsmeans_x (drop = df tvalue probt) 
			diffs = diffs_x (keep = estimate probt &iv. _&iv.);
	run;

	data diffs_&dv. (rename = (estimate = Compr_mean));
		set diffs_x;
		mergeme = 1;
		var = left(trim(put(&iv,5.0)));
		var2 = left(trim(put(_&iv,5.0)));
	run;

	data lsmeans_&dv. (rename = (estimate = lsmean effect = iv));
		length dv research_q delivery_ty $50.;
		set lsmeans_x;
		mod_number = &mod.;
		dv = "&dv.";
		mod_type = "&mod_ty.";
		research_q = "&research_q.";
		lsmean_level = left(trim(put(&iv,5.0)));
		mergeme = 1;
		%if "&ds" = "sheet3" %then %do;
        delivery_ty = "&delivery."; 
		%end;
		%else %do;
		delivery_ty = "NA Sheet1,2"; 
        %end;
	run;

	data merge_&mod. (drop=mergeme);
		length iv $50.;
		merge diffs_&dv. (drop =  &iv _&iv.) 
			lsmeans_&dv. (drop =  &iv);
		by mergeme;
	run;

	data allaov;
		retain dv mod_type mod_number research_q iv lsmean_level lsmean stderr;
		set allaov merge_&mod.;
	run;

%mend;

data allaov; set _null_; data lsmeans_x; set _null_; data diffs_x; set _null_; run;

%macro loop_aov(ds_num,iv,delivery);
%let dv1 =d_post neut plt;
%let dv2 =cmax auclast tlast;
%let dv3 =cma auc;

	%let vars = &&dv&ds_num.;
	%put vars: &vars.;
	%let word_cnt=%sysfunc(countw(&vars));

	%do i = 1 %to &word_cnt.;
		%aov(sheet&ds_num.,aov,&i,&i,&iv.,%scan(&vars.,&i),&delivery.);
	%end;
%mend;

%loop_aov(1,delivery,);
%loop_aov(2,cycle,);
%loop_aov(3,grp,emb);
%loop_aov(3,grp,ia);

sheet2
cycle iv rm- id

sheet1
delivery iv rm -id

sheet3
grp iv where delivery emb ia rm - id



%macro means_desc(ds);
	%let vars = &cont_describe.;
	%let word_cnt=%sysfunc(countw(&vars));

	%do i = 1 %to &word_cnt.;

		proc means data =  &ds. noprint;
			var %scan(&vars.,&i);
			output out= nmiss nmiss(%scan(&vars.,&i))=;
			output out= mean mean(%scan(&vars.,&i))=;
			output out= median median(%scan(&vars.,&i))=;
			output out= std std(%scan(&vars.,&i))=;
			output out= min min(%scan(&vars.,&i))=;
			output out= max max(%scan(&vars.,&i))=;
			output out= p25 p25(%scan(&vars.,&i))=;
			output out= p75 p75(%scan(&vars.,&i))=;
		run;

		data mean (rename = (%scan(&vars.,&i) = mean _freq_ = n));
			set mean (drop =  _type_);

		data nmiss (rename = %scan(&vars.,&i) = nmiss);
			set nmiss (keep = %scan(&vars.,&i));

		data median (rename = %scan(&vars.,&i) = median);
			set median (keep = %scan(&vars.,&i));

		data std (rename = %scan(&vars.,&i) = std);
			set std (keep = %scan(&vars.,&i));

		data min (rename = %scan(&vars.,&i) = min);
			set min (keep = %scan(&vars.,&i));

		data max (rename = %scan(&vars.,&i) = max);
			set max (keep = %scan(&vars.,&i));

		data p75 (rename = %scan(&vars.,&i) = p75);
			set p75 (keep = %scan(&vars.,&i));

		data p25 (rename = %scan(&vars.,&i) = p25);
			set p25 (keep = %scan(&vars.,&i));

		data all_mean_%scan(&vars.,&i);
			length ds $30.;
			merge nmiss mean median std min max p25 p75;
			ds = "%scan(&vars.,&i)";
		run;

		data all_mean;
			set all_mean all_mean_%scan(&vars.,&i);
		run;

	%end;
%mend;

data all_mean;
	set _null_;
run;

%means_desc(sunlight);

%macro freqs1(ds);
	%let vars = &cat_describe;
	%let word_cnt=%sysfunc(countw(&vars));

	%do i = 1 %to &word_cnt.;

		proc freq data =  &ds.;
			tables %scan(&vars.,&i) / chisq;
			ods output onewayfreqs = owf OneWayChiSq= owc;
		run;

		data freq (rename = f_%scan(&vars.,&i) = var_mid);
			set owf (drop = %scan(&vars.,&i));
		run;

		data chi (drop = name1 label1);
			set owc;
			where name1 = "P_PCHI";
		run;

		data freq_%scan(&vars.,&i);
			length table var_level cValue1 $50.;
			merge freq chi;
			by table;
			p_value = nvalue1;
			var_level = left(var_mid);

		data freq_tot (drop= var_mid nvalue1 cvalue1 );
			set freq_tot freq_%scan(&vars.,&i);
		run;

	%end;
%mend;

data freq_tot;
	set _null_;
run;

%freqs1(sunlight);


%let cont_describe = 
	age
	wt
	pc_pre
	hct_pre
	mcv_pre
	alb_pre
	bun_pre
;

%macro normality_test(dv=,iv=);
	ods tagsets.sasreport13(id=egsr) gtitle gfootnote;

	proc reg data=sunlight;
		model &dv. = &iv.;
		output out=residuals(keep=&iv. &dv. single r lev cd dffit)
			rstudent=r h=lev cookd=cd dffits=dffit;
	run;

	proc univariate data=residuals plot plotsize=30 normal;
		var r;
		output out=test_results_&dv probm=norm_test_&dv;
	run;

	data test_results_&dv;
		set test_results_&dv;
		variable = "&iv.          ";
	run;

	data normality_test;
		set normality_test test_results_&dv;
	run;

%mend;

data normality_test;
	set _null_;
run;

%macro loop;
	%let vars = &cont_describe;
	%let word_cnt=%sysfunc(countw(&vars));

	%do i = 1 %to &word_cnt.;
		%normality_test(dv=ih_eh,iv=%scan(&vars.,&i));
	%end;
%mend;

%loop;




%macro freq_cust(analysis,r_q_n,research,iv1,iv2);

	proc freq data =  sunlight;
		tables &iv1.*&iv2. / chisq;
		ods output CrossTabFreqs = freqs ChiSq= chi FishersExact = fish;
	run;

	run;

	data freqs_in;
		set freqs (drop =  _table_ rowpercent colpercent);
		where _type_ in ("11" "00");
		drop _type_;
	run;

	data ctf_&research. (drop =  &iv1. &iv2.);
		retain r_q_n analysis research_q;
		length r_q_n research_q $50.;
		set freqs_in;
		analysis = "&analysis.";
		research_q = "&research.";
		r_q_n = "&r_q_n.";
		iv1 = put(&iv1.,5.0);
		iv2 = put(&iv2.,5.0);
	run;

	data ctf_all;
		set ctf_all ctf_&research.;
	run;

	%if %sysfunc(exist(fish)) %then
		%do;

			data f_&research. (rename=nvalue1=P_FISH);
				set fish (keep =  table name1 nvalue1);
				where name1 in ("P_TABLE");
				drop name1;
			run;

			data f_p;
				set f_p f_&research.;
			run;

		%end;

	%if %sysfunc(exist(chi)) %then
		%do; 
			%if (&iv1. = Br OR &iv1. = size OR &iv1. = Hx_pneum OR &iv1. = Skel) OR (&iv2. = abx_cor1) %then %do;
				data chi_&research. (rename=prob=P_CHI);
					set chi (keep = table statistic prob);
					where statistic in ("Mantel-Haenszel Chi-Square");
					Research_Question = "&research";
				run;

				data chi_p;
					set chi_p chi_&research.;
				run;
			%end;
			%else %do;
				data chi_&research. (rename=prob=P_CHI);
					set chi (keep =  table statistic prob);
					where statistic in ("Chi-Square");
					length research_question $10;
					Research_Question = "&research";
				run;

				data chi_p;
					set chi_p chi_&research.;
				run;

			%end;

		%end;

	proc datasets library=work nolist nodetails;
		delete chi fish freq p_&research. chi_&research.;
	run;

%mend;

data ctf_all;
	set _null_;

data f_p;
	set _null_;

data chi_p;
	set _null_;
run;

%let cat_describe = 
	size 
	male
;

%freq_cust(sunlight,1,1,size,ih_eh);
%freq_cust(sunlight,2,2,male,ih_eh);



%let cont_describe = 
pcf_post
pcf_tot
pcf_per
hcf_post
hcf_tot
hcf_per
mcf_post
mcf_tot
mcf_per
alf_post
alf_tot
alf_per
buf_post
buf_tot
buf_per
;

%macro means_desc(ds);
	%let vars = &cont_describe.;
	%let word_cnt=%sysfunc(countw(&vars));

	%do i = 1 %to &word_cnt.;

		proc means data =  &ds. noprint;
			var %scan(&vars.,&i);
			output out= nmiss nmiss(%scan(&vars.,&i))=;
			output out= mean mean(%scan(&vars.,&i))=;
			output out= median median(%scan(&vars.,&i))=;
			output out= std std(%scan(&vars.,&i))=;
			output out= min min(%scan(&vars.,&i))=;
			output out= max max(%scan(&vars.,&i))=;
			output out= p25 p25(%scan(&vars.,&i))=;
			output out= p75 p75(%scan(&vars.,&i))=;
		run;

		data mean (rename = (%scan(&vars.,&i) = mean _freq_ = n));
			set mean (drop =  _type_);

		data nmiss (rename = %scan(&vars.,&i) = nmiss);
			set nmiss (keep = %scan(&vars.,&i));

		data median (rename = %scan(&vars.,&i) = median);
			set median (keep = %scan(&vars.,&i));

		data std (rename = %scan(&vars.,&i) = std);
			set std (keep = %scan(&vars.,&i));

		data min (rename = %scan(&vars.,&i) = min);
			set min (keep = %scan(&vars.,&i));

		data max (rename = %scan(&vars.,&i) = max);
			set max (keep = %scan(&vars.,&i));

		data p75 (rename = %scan(&vars.,&i) = p75);
			set p75 (keep = %scan(&vars.,&i));

		data p25 (rename = %scan(&vars.,&i) = p25);
			set p25 (keep = %scan(&vars.,&i));

		data all_mean_%scan(&vars.,&i);
			length ds $30.;
			merge nmiss mean median std min max p25 p75;
			ds = "%scan(&vars.,&i)";
		run;

		data all_mean;
			set all_mean all_mean_%scan(&vars.,&i);
		run;

	%end;
%mend;

data all_mean;
	set _null_;
run;

%means_desc(sunlight(where=(ih_eh=1)));
%means_desc(sunlight(where=(ih_eh=2)));


***** ANOVA Analysis;
%macro aov(mod_ty,mod,research_q,iv,dv);
	/*ods output lsmeans tests3 diffs;*/
	proc mixed data=sunlight;
		class &iv.;
		model &dv. = &iv.  / solution;
		lsmeans &iv.  / pdiff=all;
		ods output lsmeans=lsmeans_x (drop = df tvalue probt) 
			diffs = diffs_x (keep = estimate probt &iv. _&iv.);
	run;

	data diffs_&dv. (rename = (estimate = Compr_mean));
		set diffs_x;
		mergeme = 1;
		var = left(trim(put(&iv,5.0)));
		var2 = left(trim(put(_&iv,5.0)));
	run;

	data lsmeans_&dv. (rename = (estimate = lsmean effect = iv));
		length dv research_q $50.;
		set lsmeans_x;
		mod_number = &mod.;
		dv = "&dv.";
		mod_type = "&mod_ty.";
		research_q = "&research_q.";
		lsmean_level = left(trim(put(&iv,5.0)));
		mergeme = 1;
	run;

	data merge_&mod. (drop=mergeme);
		length iv $50.;
		merge diffs_&dv. (drop =  &iv _&iv.) 
			lsmeans_&dv. (drop =  &iv);
		by mergeme;
	run;

	data allaov;
		retain dv mod_type mod_number research_q iv lsmean_level lsmean stderr;
		set allaov merge_&mod.;
	run;

%mend;

data allaov;
	set _null_;

data lsmeans_x;
	set _null_;

data diffs_x;
	set _null_;
run;

%macro loop_aov;
	%let vars = &cont_describe;
	%let word_cnt=%sysfunc(countw(&vars));

	%do i = 1 %to &word_cnt.;
		%aov(aov,&i,&i,ih_eh,%scan(&vars.,&i));
	%end;
%mend;

%loop_aov;

***** ANOVA Analysis;
%macro aov_subset(mod_ty,mod,research_q,iv,dv);
	/*ods output lsmeans tests3 diffs;*/
	proc mixed data=sunlight(where=(ih_eh=1));
		class &iv.;
		model &dv. = &iv.  / solution;
		lsmeans &iv.  / pdiff=all;
		ods output lsmeans=lsmeans_x (drop = df tvalue probt) 
			diffs = diffs_x (keep = estimate probt &iv. _&iv.);
	run;

	data diffs_&dv. (rename = (estimate = Compr_mean));
		set diffs_x;
		mergeme = 1;
		var = left(trim(put(&iv,5.0)));
		var2 = left(trim(put(_&iv,5.0)));
	run;

	data lsmeans_&dv. (rename = (estimate = lsmean effect = iv));
		length dv research_q $50.;
		set lsmeans_x;
		mod_number = &mod.;
		dv = "&dv.";
		mod_type = "&mod_ty.";
		research_q = "&research_q.";
		lsmean_level = left(trim(put(&iv,5.0)));
		mergeme = 1;
	run;

	data merge_&mod. (drop=mergeme);
		length iv $50.;
		merge diffs_&dv. (drop =  &iv _&iv.) 
			lsmeans_&dv. (drop =  &iv);
		by mergeme;
	run;

	data allaov;
		retain dv mod_type mod_number research_q iv lsmean_level lsmean stderr;
		set allaov merge_&mod.;
	run;

%mend;

%macro loop_aov;
	%let vars = &cont_describe;
	%let word_cnt=%sysfunc(countw(&vars));

	%do i = 1 %to &word_cnt.;
		%aov_subset(aov,&i,&i,ih_eh,%scan(&vars.,&i));
	%end;
%mend;

%loop_aov;

%macro freq_cust(analysis,r_q_n,research,iv1,iv2);

	proc freq data =  sunlight;
		tables &iv1.*&iv2. / chisq;
		ods output CrossTabFreqs = freqs ChiSq= chi FishersExact = fish;
	run;

	run;

	data freqs_in;
		set freqs (drop =  _table_ rowpercent colpercent);
		where _type_ in ("11" "00");
		drop _type_;
	run;

	data ctf_&research. (drop =  &iv1. &iv2.);
		retain r_q_n analysis research_q;
		length r_q_n research_q $50.;
		set freqs_in;
		analysis = "&analysis.";
		research_q = "&research.";
		r_q_n = "&r_q_n.";
		iv1 = put(&iv1.,5.0);
		iv2 = put(&iv2.,5.0);
	run;

	data ctf_all;
		set ctf_all ctf_&research.;
	run;

	%if %sysfunc(exist(fish)) %then
		%do;

			data f_&research. (rename=nvalue1=P_FISH);
				set fish (keep =  table name1 nvalue1);
				where name1 in ("P_TABLE");
				drop name1;
			run;

			data f_p;
				set f_p f_&research.;
			run;

		%end;

	%if %sysfunc(exist(chi)) %then
		%do; 
			%if (&iv1. = Br OR &iv1. = size OR &iv1. = Hx_pneum OR &iv1. = Skel) OR (&iv2. = abx_cor1) %then %do;
				data chi_&research. (rename=prob=P_CHI);
					set chi (keep = table statistic prob);
					where statistic in ("Mantel-Haenszel Chi-Square");
					Research_Question = "&research";
				run;

				data chi_p;
					set chi_p chi_&research.;
				run;
			%end;
			%else %do;
				data chi_&research. (rename=prob=P_CHI);
					set chi (keep =  table statistic prob);
					where statistic in ("Chi-Square");
					length research_question $10;
					Research_Question = "&research";
				run;

				data chi_p;
					set chi_p chi_&research.;
				run;

			%end;

		%end;

	proc datasets library=work nolist nodetails;
		delete chi fish freq p_&research. chi_&research.;
	run;

%mend;

data ctf_all;
	set _null_;

data f_p;
	set _null_;

data chi_p;
	set _null_;
run;

%let cat_describe = 
	size 
	male
;

%freq_cust(sunlight,1,1,size,ih_eh);
%freq_cust(sunlight,2,2,male,ih_eh);


%aov_subset(aov,10,10,oc,pc_pre);
%aov_subset(aov,10,10,oc,hct_pre);
%aov_subset(aov,10,10,oc,mcv_pre);
%aov_subset(aov,10,10,oc,alb_pre);
%aov_subset(aov,10,10,oc,bun_pre);

%aov_subset(aov,11,11,oc,pcf_tot);
%aov_subset(aov,11,11,oc,pcf_per);
%aov_subset(aov,11,11,oc,hcf_tot);
%aov_subset(aov,11,11,oc,hcf_per);
%aov_subset(aov,11,11,oc,mcf_tot);
%aov_subset(aov,11,11,oc,mcf_per);
%aov_subset(aov,11,11,oc,alf_tot);
%aov_subset(aov,11,11,oc,alf_per);
%aov_subset(aov,11,11,oc,buf_tot);
%aov_subset(aov,11,11,oc,buf_per);

