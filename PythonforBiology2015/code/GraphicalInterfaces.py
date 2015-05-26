import re
try:
  import tkinter
  from tkinter import filedialog, messagebox
except:
  import Tkinter as tkinter
  import tkFileDialog as filedialog
  import tkMessageBox as messagebox

from Sequences import proteinTranslation, STANDARD_GENETIC_CODE

try:
  from Bio import SeqIO
  
except ImportError:
  print('\n* * * * * BioPython not installed * * * * * * * * ')
  print('\n* * * * * Some functions will not work  * * * * * ')


from PySide import QtCore, QtGui 
# or from PyQt4 import QtCore, QtGui 

rootWindow = tkinter.Tk()

label = tkinter.Label(rootWindow, text='Hello World')

label.pack()

rootWindow.mainloop()


class SequenceTkGui(tkinter.Tk):

  def __init__(self):

    tkinter.Tk.__init__(self)

    self.grid_columnconfigure(5, weight=1)
    self.grid_rowconfigure(1, weight=1)
    self.grid_rowconfigure(4, weight=1)
 
    self.label1 = tkinter.Label(self, text='Enter 1-Letter DNA Sequence:')
    self.label1.grid(row=0, column=0, columnspan=6, sticky=tkinter.EW)

    self.seqTextBox = tkinter.Text(self)
    self.seqTextBox.grid(row=1, column=0, columnspan=6, 
                         sticky=tkinter.NSEW)

    self.clearButton = tkinter.Button(self, text='Clear',
                                      command=self.clearSeq)
    self.clearButton.grid(row=2, column=0, sticky=tkinter.W)

    self.loadButton = tkinter.Button(self, text='Load FASTA',
                                     command=self.loadFasta)
    self.loadButton.grid(row=2, column=1, sticky=tkinter.W)

    self.transButton = tkinter.Button(self, text='Translate',
                                      command=self.seqTranslate)
    self.transButton.grid(row=2, column=2, sticky=tkinter.W)

    self.compButton = tkinter.Button(self, text='Composition',
                                     command=self.seqComposition)
    self.compButton.grid(row=2, column=3, sticky=tkinter.W)

    self.findButton = tkinter.Button(self, text='Find:',
                                     command=self.seqFind)
    self.findButton.grid(row=2, column=4, sticky=tkinter.EW)

    self.findEntry = tkinter.Entry(self)
    self.findEntry.grid(row=2, column=5, sticky=tkinter.EW)

    self.label2 = tkinter.Label(self, text='Text output:')
    self.label2.grid(row=3, column=0, columnspan=6, sticky=tkinter.W)

    self.outTextBox = tkinter.Text(self)
    self.outTextBox.grid(row=4, column=0, columnspan=6,             
                         sticky=tkinter.NSEW)

    self.closeButton = tkinter.Button(self, text='Quit',                          
                                      command=self.destroy)
    self.closeButton.grid(row=5, column=5, sticky=tkinter.EW)
    self.closeButton.config(bg='yellow')

  def clearSeq(self):

    self.seqTextBox.delete('0.0', tkinter.END)

  def setSequence(self, text):

    self.clearSeq()
    self.seqTextBox.insert(tkinter.END, text)
 
  def getSequence(self):

    seq = self.seqTextBox.get('0.0', tkinter.END)
    seq = re.sub('\s+','',seq)
    seq = seq.upper()

    return seq
   
  def showText(self, text):

    if text[-1] != '\n':
      text += '\n'
    self.outTextBox.insert(tkinter.END, text)

  def clearOutput(self):

    self.outTextBox.delete('0.0', tkinter.END)

  def loadFasta(self):
 
    fileObj = filedialog.askopenfile(parent=self, mode='rU', 
                                     title='Choose a FASTA file')

    if fileObj:
      from Bio import SeqIO
      for entry in SeqIO.parse(fileObj, 'fasta'):
        self.setSequence(entry.seq)
        break

      fileObj.close()

  def seqTranslate(self):

    seq = self.getSequence()

    self.clearOutput()
    self.showText('DNA sequence')
    self.showText(seq)
    self.showText('Protein sequence')
 
    for indent in range(3):

      proteinSeq = proteinTranslation(seq[indent:], STANDARD_GENETIC_CODE)
      proteinSeq = ''.join(proteinSeq)
      spaces = ' ' * indent
 
      text = 'Reading frame %d\n%s%s' % (indent, spaces, proteinSeq)

      self.showText(text)

  def seqComposition(self):

    self.clearOutput()
    seq = self.getSequence()

    n = 0.0
    counts = {}
    for letter in seq:
      counts[letter] = counts.get(letter, 0) + 1
      n += 1.0
 
    letters = counts.keys()
    letters.sort()
 
    text = "Composition:"
    for letter in letters:
      text += ' %s;%.2f%%' % (letter, counts[letter] * 100 / n)
 
    self.showText(text)

  def seqFind(self):
 
    self.clearOutput()
    
    query = self.findEntry.get()
    query = query.strip()
    
    if not query:
      messagebox.showwarning("Warning", "Search sequence was blank")
      return
      
    seq  = self.getSequence()

    if query in seq:
      text = "Locations of %s" % (query)
      self.showText(text)
      win = len(query)
      
      for i in range(len(seq)-win):
        if seq[i:i+win] == query:
          self.showText(' %d' % i)
    
    else:
      text = "Sub-sequence %s not found" % (query)
      self.showText(text)


if __name__ == '__main__':
 
  window = SequenceTkGui()
  window.mainloop()


class SequenceQtGui(QtGui.QWidget):

  def __init__(self):

    QtGui.QWidget.__init__(self, parent=None)
    
    grid = QtGui.QGridLayout(self)
    
    grid.setColumnStretch(5, 1)
    grid.setRowStretch(1, 1)
    grid.setRowStretch(4, 1)
    
    leftRight = QtCore.Qt.AlignLeft | QtCore.Qt.AlignRight

    self.label1 = QtGui.QLabel(text='Enter 1-Letter DNA Sequence:',                          
                               parent=self)
    grid.addWidget(self.label1, 0, 0)

    self.seqTextBox = QtGui.QPlainTextEdit(parent=self)
    grid.addWidget(self.seqTextBox, 1, 0, 1, 6, align=leftRight)

    self.clearButton = QtGui.QPushButton(text='Clear', parent=self)
    self.clearButton.clicked.connect(self.clearSeq)
    grid.addWidget(self.clearButton, 2, 0)

    #PyQt4 uses:
    #self.connect(self, QtCore.SIGNAL('clicked()'), self.clearSeq)

    self.loadButton = QtGui.QPushButton(text='Load FASTA', parent=self)
    self.loadButton.clicked.connect(self.loadFasta)
    grid.addWidget(self.loadButton, 2, 1)

    self.transButton = QtGui.QPushButton(text='Translate', parent=self)
    self.transButton.clicked.connect(self.seqTranslate)
    grid.addWidget(self.transButton, 2, 2)

    self.compButton = QtGui.QPushButton(text='Composition', parent=self)
    self.compButton.clicked.connect(self.seqComposition)
    grid.addWidget(self.compButton, 2, 3)

    self.findButton = QtGui.QPushButton(text='Find:', parent=self)
    self.findButton.clicked.connect(self.seqFind)
    grid.addWidget(self.findButton, 2, 4)

    self.findEntry = QtGui.QLineEdit(parent=self)
    grid.addWidget(self.findEntry, 2, 5)

    self.label2 = QtGui.QLabel(text='Text output:', parent=self)
    grid.addWidget(self.label2, 3, 0, 1, 6)

    self.outTextBox = QtGui.QPlainTextEdit(parent=self)
    self.outTextBox.setLineWrapMode(QtGui.QPlainTextEdit.NoWrap)
    grid.addWidget(self.outTextBox, 4, 0, 1, 6, align=leftRight)

    self.closeButton = QtGui.QPushButton(self, text='Quit')
    self.closeButton.clicked.connect(self.destroy)
    grid.addWidget(self.closeButton, 5, 5, align=leftRight)

  def clearSeq(self):

    self.seqTextBox.clear()

  def setSequence(self, text):

    self.seqTextBox.setPlainText(text)
    
  def showText(self, text):

    self.outTextBox.appendPlainText(text)

  def clearOutput(self):

    self.outTextBox.clear()

  def getSequence(self):

    seq = self.seqTextBox.toPlainText()
    seq = re.sub('\s+','',seq)
    seq = seq.upper()

    return seq

  def loadFasta(self):
    
    msg = 'Choose a FASTA file'
    filePath, filtr = QtGui.QFileDialog.getOpenFileName(self, msg)
    
    if filePath: # Something was selected
      fileObj = open(filePath, 'rU')

      from Bio import SeqIO
      for entry in SeqIO.parse(fileObj, 'fasta'):
        self.setSequence(str(entry.seq))
        break

      fileObj.close()

  def seqTranslate(self):

    seq = self.getSequence()

    self.clearOutput()
    self.showText('DNA sequence')
    self.showText(seq)
    self.showText('Protein sequence')
 
    for indent in range(3):

      proteinSeq = proteinTranslation(seq[indent:], STANDARD_GENETIC_CODE)
      proteinSeq = ''.join(proteinSeq)
      spaces = ' ' * indent
 
      text = 'Reading frame %d\n%s%s' % (indent, spaces, proteinSeq)

      self.showText(text)

  def seqComposition(self):

    self.clearOutput()
    seq = self.getSequence()

    n = 0.0
    counts = {}
    for letter in seq:
      counts[letter] = counts.get(letter, 0) + 1
      n += 1.0
 
    letters = counts.keys()
    letters.sort()
 
    text = "Composition:"
    for letter in letters:
      text += ' %s;%.2f%%' % (letter, counts[letter] * 100 / n)
 
    self.showText(text)

  def seqFind(self):
 
    self.clearOutput()
    
    query = self.findEntry.text()
    query = query.strip()
    
    if not query:
      QtGui.QMessageBox.warning(self, "Warning",
                                "Search sequence was blank")
      return
      
    seq = self.getSequence()
    
    self.seqTextBox.find(query)
    
    if query in seq:
      text = "Locations of %s" % (query)
      self.showText(text)
      win = len(query)
      
      for i in range(len(seq)-win):
        if seq[i:i+win] == query:
          self.showText(' %d' % i)
    
    else:
      text = "Sub-sequence %s not found" % (query)
      self.showText(text)
  
if __name__ == '__main__':
  
  import sys
  
  app = QtGui.QApplication(['Qt Sequence Example'])
  
  window = SequenceQtGui()
  window.show()
  
  sys.exit(app.exec_())

