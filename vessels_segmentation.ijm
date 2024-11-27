setBatchMode(false);

run("CLIJ2 Macro Extensions", "cl_device=[NVIDIA RTX A5000]");
Ext.CLIJ2_clear();

// 1. MEDIAN FILTER
image14 = Image.title;
Ext.CLIJ2_push(image14);
image15 = "median";
radiusx = 2.0;
radiusy = 2.0;
radiusz = 1.0;
Ext.CLIJ2_median3DSphere(image14, image15, radiusx, radiusy, radiusz);
Ext.CLIJ2_pull(image15);

// 2. DIFFERENCE OF GAUSSIAN FILTER
Ext.CLIJ2_push(image15);
image16 = "dog";
sigma1x = 3.0; // To optimize
sigma1y = 3.0; // Same than sigma1x
sigma2x = 6.0; // To optimize
sigma2y = 6.0; // Same than sigma2x
Ext.CLIJ2_differenceOfGaussian2D(image15, image16, sigma1x, sigma1y, sigma2x, sigma2y);
Ext.CLIJ2_pull(image16);

// 3. THRESHOLD
Ext.CLIJ2_push(image16);
image17 = "threshold";
Ext.CLIJ2_thresholdTriangle(image16, image17); // To optimize
Ext.CLIJ2_pull(image17);

// 4. CLOSING FILTER
radiusx = 4.0;
radiusy = 4.0;
radiusz = 1.0;

// 4.a. Maximum filter
Ext.CLIJ2_push(image17);
image18 = "maximum";
Ext.CLIJ2_maximum3DSphere(image17, image18, radiusx, radiusy, radiusz);
Ext.CLIJ2_pull(image18);

// 4.b. Minimum filter
Ext.CLIJ2_push(image18);
image19 = "minimum";
Ext.CLIJ2_minimum3DSphere(image18, image19, radiusx, radiusy, radiusz);
Ext.CLIJ2_pull(image19);

// 5. MEDIAN FILTER
Ext.CLIJ2_push(image19);
image20 = "result";
radiusx = 2.0;
radiusy = 2.0;
radiusz = 1.0;
Ext.CLIJ2_median3DSphere(image19, image20, radiusx, radiusy, radiusz);
Ext.CLIJ2_pull(image20);

setThreshold(1, 255);
run("Make Binary", "black");

close("median");
close("dog");
close("threshold");
close("maximum");
close("minimum");
Ext.CLIJ2_clear();

setBatchMode(false);
