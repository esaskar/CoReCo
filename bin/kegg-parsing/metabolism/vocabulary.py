# volabulary.py - metabolism.vocabulary
#
#

class Vocabulary:

    def __init__(self):
        self.ids = {}
        self.synonyms = {}
        self.commonNames = {}

    def addSynonym(self, name, synonym, commonName = 0):
        self.ids[name] = 1

        if commonName:
            self.commonNames[synonym] = 1
        
        if name not in self.synonyms:
            self.synonyms[name] = {}

        if synonym not in self.synonyms:
            self.synonyms[synonym] = {}

        if synonym not in self.synonyms[name]:
            self.synonyms[name][synonym] = 1

        if name not in self.synonyms[synonym]:
            self.synonyms[synonym][name] = 1

    def getSynonyms(self, name):
        if name in self.synonyms:
            return self.synonyms[name].keys()
        else:
            return []

    def getCommonName(self, name):
        if name in self.synonyms:
            for k in self.synonyms[name]:
                if k in self.commonNames:
                    return k
        return None

    def getNonIdSynonym(self, name):
        if name in self.synonyms:
            keys = self.synonyms[name].keys()
            for k in keys:
                if k not in self.ids:
                    return k
        return None
                
    def getIdForSynonym(self, synonym):
        if synonym in self.ids:
            return synonym

        if synonym not in self.synonyms:
            return None
        
        for name in self.synonyms[synonym]:
            if name in self.ids:
                return name

        return None

    def __str__(self):
        return "Vocabulary with %d keys" % (len(self.synonyms))
