lafs = javax.swing.UIManager.getInstalledLookAndFeels
for lafIdx = 1:length(lafs),  disp(lafs(lafIdx));  end
javax.swing.UIManager.setLookAndFeel('com.sun.java.swing.plaf.windows.WindowsLookAndFeel');
