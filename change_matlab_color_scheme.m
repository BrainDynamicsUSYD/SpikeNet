% change MATLAB color scheme
% from http://www.mikesoltys.com/2013/02/08/matlab-tip-change-the-color-scheme-to-be-easier-on-your-eyes/

% Also see
% http://ethanschoonover.com/solarized


com.mathworks.services.Prefs.setBooleanPref('ColorsUseSystem',0);
com.mathworks.services.Prefs.setColorPref('ColorsBackground',java.awt.Color(7/255,54/255,66/255));
com.mathworks.services.Prefs.setColorPref('ColorsText',java.awt.Color(133/255,153/255,0/255));
com.mathworks.services.Prefs.setColorPref('Colors_M_Comments',java.awt.Color(181/255,137/255,0/255));
com.mathworks.services.Prefs.setColorPref('Colors_M_Keywords',java.awt.Color(38/255,139/255,210/255));
com.mathworks.services.Prefs.setColorPref('Colors_M_Strings',java.awt.Color(42/255,161/255,152/255));