# -*- coding: utf-8 -*-
"""
Created on Fri Sep 26 16:28:14 2014

@author: Nicky Jay feat Angus M
Welcome to Clicky 2.0"""

import matplotlib.pyplot as plt
from matplotlib.widgets import Cursor, Button
plt.rcParams['font.size']=8
from PIL import Image
import os
import mayavi.mlab as mlab
from ConfigParser import SafeConfigParser
import numpy as np
import scipy.ndimage as ndi
import sys

#Type in z axis slice locations


class viewer_2d(object):
    
    def __init__(self,top,mid1,mid2,mid3,bottom):
        """
        Shows a given array in a 2d-viewer.
        Input: 2d arrays. For each slice level selected
        Also where the buttons on the viewer are added
        """
        
        self.top = top
        self.fig=plt.figure(figsize=(12,8))
        #Doing some layout with subplots:
        self.fig.subplots_adjust(0.05,0.05,0.98,0.98,0.1)
        self.top_im = plt.subplot2grid((7,2),(0,0),rowspan=1,colspan=2)
        self.top_im.imshow(self.top)
        self.top_im.autoscale(1,'both',1)
        
        self.mid1 = mid1
        self.mid1_im=plt.subplot2grid((7,2),(1,0),rowspan=1,colspan=2)
        self.mid1_im.imshow(self.mid1)
        self.mid1_im.autoscale(1,'both',1)
        
        self.mid2 = mid2
        self.mid2_im=plt.subplot2grid((7,2),(2,0),rowspan=1,colspan=2)
        self.mid2_im.imshow(self.mid2)
        self.mid2_im.autoscale(1,'both',1)
        
        self.mid3 = mid3
        self.mid3_im=plt.subplot2grid((7,2),(3,0),rowspan=1,colspan=2)
        self.mid3_im.imshow(self.mid3)
        self.mid3_im.autoscale(1,'both',1)
        
        self.bottom = bottom
        self.bottom_im=plt.subplot2grid((7,2),(4,0),rowspan=1,colspan=2)
        self.bottom_im.imshow(self.bottom)
        self.bottom_im.autoscale(1,'both',1)
        
        # internal storage for coordinates
        self.top_list    = []
        self.mid1_list    = []
        self.mid2_list    = []
        self.mid3_list    = []
        self.bottom_list = []
        self.last = None
        
        self.cell_number = 1
        #Adding widgets, to not be gc'ed, they are put in a list:        
        top_cursor = Cursor(self.top_im, useblit=True, color='black', linewidth=2 )
        mid1_cursor = Cursor(self.mid1_im, useblit=True, color='black', linewidth=2 )
        mid2_cursor = Cursor(self.mid2_im, useblit=True, color='black', linewidth=2 )
        mid3_cursor = Cursor(self.mid3_im, useblit=True, color='black', linewidth=2 )
        bottom_cursor = Cursor(self.bottom_im, useblit=True, color='black', linewidth=2 )
        but_ax = plt.subplot2grid((7,3),(5,0),colspan=1)
        undo_button = Button(but_ax,'Undo')
        but_ax2 = plt.subplot2grid((7,3),(5,2),colspan=1)
        save_button = Button(but_ax2,'Save')
        but_ax3 = plt.subplot2grid((7,3),(5,1),colspan=1)
        next_button = Button(but_ax3, 'Next Cell')
        self._widgets = [top_cursor, bottom_cursor, mid1_cursor,mid2_cursor,mid3_cursor, undo_button, save_button, next_button]
        
        # connect events
        undo_button.on_clicked(self.delete_last)
        save_button.on_clicked(self.save)
        next_button.on_clicked(self.next_cell)
        self.fig.canvas.mpl_connect('button_press_event',self.click)
 
    def next_cell(self, evt):
        self.cell_number += 1
        self.top_list    = []
        self.mid1_list    = []
        self.mid2_list    = []
        self.mid3_list    = []
        self.bottom_list = []
        print "Next cell selected"
        
    def delete_last(self, evt):
        del self.last[-1]
      
    def save(self, evt):
        num_cells = len(self.top_list)
        parser = SafeConfigParser()        
        parser.read(r'C:\Users\Dell\Desktop\curved_config\P6_5_ribbon.ini')
        
        """if section Cell 1 does not exist, make it"""
        if not parser.has_section('Cell %s' % self.cell_number):
            parser.add_section('Cell %s' % self.cell_number)
        if not parser.has_section('misc'):
                parser.add_section('misc')
            
        parser.set('misc', 'number_of_clicks', str(num_cells))            
        parser.set('Cell %s' % self.cell_number, 'Top',    str((self.top_list) + (self.mid1_list) +(self.mid2_list)+(self.mid3_list) + (self.bottom_list)))
        
        out = open(r'C:\Users\Dell\Desktop\curved_config\P6_5_ribbon.ini', 'w')
        parser.write(out)
        print "Data saved to config file"
      
      
    def click(self, event):
        """
        What to do, if a click on the figure happens:
            1. Check which axis
            2. Get data coord's.
            3. Plot resulting data.
            4. Update Figure
        """
        
        if event.inaxes in [self.top_im, self.mid1_im, self.mid2_im, self.mid3_im, self.bottom_im]:
            if event.inaxes == self.top_im:
                self.last = self.top_list
                self.y_pos = top_y
            if event.inaxes == self.mid1_im:
                self.last = self.mid1_list
                self.y_pos = mid1_y
            if event.inaxes == self.mid2_im:
                self.last = self.mid2_list
                self.y_pos = mid2_y
            if event.inaxes == self.mid3_im:
                self.last = self.mid3_list
                self.y_pos = mid3_y
            if event.inaxes == self.bottom_im:
                self.last = self.bottom_list
                self.y_pos = bottom_y
            
          
            # Get nearest data
            xpos = event.xdata
            zpos = event.ydata

            # Check which mouse button:
            if event.button == 1:
                print zpos, self.y_pos, xpos
                self.last.append([zpos, self.y_pos, xpos])

            elif event.button==3:
                #Plot it
                pass

        #Show it
        plt.draw()

"""Different levels for slicing through the image (in pixels) """     
top_y = 340
mid1_y = 390
mid2_y = 440
mid3_y = 480
bottom_y = 520
parser = SafeConfigParser()
"""location of confocal images""" 
wdir = os.path.normpath(r'C:\Users\DELL\Desktop\Confocal images\Glyofixx\P6')
"""location of configfiles"""
parser.read('C:\Users\DELL\Desktop\curved_config\P6_5_ribbon.ini')
mtiff_file = parser.get('section1','mask_file')
im = os.path.join(wdir, mtiff_file)
nimg = parser.getint('section1', 'number_of_slices')
img = Image.open(im)
img.seek(0)
nx, ny = img.size
big_arr = np.zeros([nimg, nx, ny])
for i in range(nimg):
    img.seek(i)
    big_arr[i] = np.array(img)

if __name__=='__main__':
  
    fig_v = viewer_2d(big_arr[:,top_y,:], big_arr[:,mid1_y,:],big_arr[:,mid2_y,:],big_arr[:,mid3_y,:], big_arr[:,bottom_y,:])
    
    #Show it
    plt.show()