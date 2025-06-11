# Grasshopper-Power-BI-Stream
![](images/Monash_University_logo.png)
## Preliminary Work

This project's functionality is aimed at generating real-time models in Power BI using the UI within Power BI. The method involves first generating data files through Power BI, then using Grasshopper to read these data files, uploading them to the Speckle website, and finally displaying them in Power BI. The plugins used in this project include Speckle and Pancake.

### Plugin Setup

Double-click manager.exe to install Speckle Manager. Then, use the Manager to download the Power BI, Rhino and Grasshopper connectors.

![](images/Connection.png)

Open Grasshopper, open Component Folder here, then drag Pancake.gha into the Component Folder shown in the figure.

![](images/GH1.png)

Open Power BI, click the three dots in the Visualizations panel, select "Import a visual from a file", and choose Speckle Power BI 3D Visual.pbiviz.

![](images/PBI1.png)

### Path Setup

After installing the plugins, open the Tunnel_Lining.3dm file in Rhino, launch Grasshopper, and load the Tunnel_Generation.gh file. In the Import Data group, Create and Bake Lining group, and Property Calculation group, update the file path to the current location of your file.

![](images/GH2.png)
![](images/GH3.png)
![](images/GH4.png)

 Open PBI UI.pbix in Power BI, click the blue-bordered area in the bottom-right corner, and update the path to the current file location.
 
![](images/PBIPY.png)

 Then, register on the Speckle website (www.speckle.systems), create a new project.
 
![](images/Speckle%20web.png)

 Create a model within this project and copy the model's URL.
 
![](images/ModelURL.png)

 Enter the URL in the following two places: Open PBI UI.pbix, click on "Transform data", and paste your URL into the indicated field in figure below. 
 
![](images/PBITD.png)

 Then click on "Data source settings" in "Transform data", click "Change source", then change the folder path of the file.
 
![](images/PBITD2.png)

 Open Tunnel_Lining.3dm, and open Tunnel_Generation.gh. Locate the "export data" group, and paste the URL into the designated field.
 
![](images/Rhino%20Web.png)

 With this, the setup is complete and you’re ready to start using the system.
 
## Steps for Use

First, open Tunnel_Lining.3dm and Tunnel_Generation.gh in the background, then launch PBI UI.pbix. Follow the steps shown below:

![](images/PBIUse.png)

 **Hint:** The current version of the Speckle plugin has a bug where repeated refreshes may cause previous models to leave behind white translucent lines. The workaround is to uncheck the URL in this visual, and then check it again.
 
![](images/PBIBUG.png)

 Then, go to the second page. On this page, you can enter parameters related to tunnel details, view various tunnel properties, and finally check the estimated carbon emissions.
 
![](images/Page%202.png)
