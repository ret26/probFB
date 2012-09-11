
function test_suite = test_plot_pSTFT
  initTestSuite;

   
function test_runs

dispFigs = 0;
lamx = [0.8;0.9];
varx = 1./(1-lamx.^2);
om = [pi/5;pi/10];

figH1 = plot_pSTFT(varx,om,lamx);

pause(0.2)

figH2 = plot_pSTFT(varx,om,lamx,1000);

pause(0.2)

figH3 = plot_pSTFT(varx,om,lamx,1,1);

pause(0.2)

if dispFigs==0
  close(figH1,figH2,figH3)
end
