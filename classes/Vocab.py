class Vocab:
    idxcnt=0
    instances={}

    @staticmethod
    def get(vtype):
        if vtype not in ('root', 'frag', 'link', 'elem'):
            raise ValueError('vtype must be \'root\', \'frag\', \'link\' or \'elem\'.')
        if vtype not in Vocab.instances:
            Vocab.instances[vtype] = Vocab(vtype)
        return Vocab.instances[vtype]
    
    @staticmethod
    def save():
        return dict(
            idxcnt=Vocab.idxcnt,
            instances=Vocab.instances
        )
    
    @staticmethod
    def load(data):
        Vocab.instances = data['instances']
        Vocab.idxcnt = data['idxcnt']
    
    @staticmethod
    def size(vtype='all'):
        if vtype=='all':
            return Vocab.idxcnt
        else:
            return Vocab.get(vtype).getSize()

    def __init__(self, vtype):
        self.type = vtype
        self.data=[]
        self.identifiers=[]
        self.indexes=[]
    
    def add(self, identifier, data=[]):
        if self.exists(identifier):
            return False
        self.identifiers.append(identifier)
        self.data.append(data)
        Vocab.idxcnt+=1
        self.indexes.append(Vocab.idxcnt)
        return True

    def exists(self, identifier):
        return identifier in self.identifiers
    
    def getIndex(self, identifier):
        if self.exists(identifier):
            return self.indexes[self.identifiers.index(identifier)]
    
    def getIdentifierByIndex(self, index):
        return self.identifiers[self.indexes.index(index)]
    
    def getDataByIndex(self, index):
        return self.data[self.indexes.index(index)]
    
    def setDataByIndex(self, index, data):
        self.data[self.indexes.index(index)]=data
        return True
    
    def getData(self, identifier):
        return self.getDataByIndex(self.getIndex(identifier))
    
    def setData(self, identifier, data):
        if not self.exists(identifier):
            return False
        else:
            return self.setDataByIndex(self.getIndex(identifier), data)
    
    def getType(self):
        return self.type

    def getSize(self):
        return len(self.indexes)