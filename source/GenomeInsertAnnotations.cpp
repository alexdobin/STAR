#include "GenomeInsertAnnotations.h"
#include "Parameters.h"
#include "ErrorWarning.h"
#include "streamFuns.h"

#include <cstdio>
#include <functional>
#include <sys/stat.h>

namespace {
const vector<string> annotationSidecarFiles = {
    "sjdbInfo.txt",
    "sjdbList.out.tab",
    "sjdbList.fromGTF.out.tab",
    "exonInfo.tab",
    "exonGeTrInfo.tab",
    "geneInfo.tab",
    "transcriptInfo.tab"
};

bool fileExists(const string &fileName)
{
    struct stat statBuf;
    return stat(fileName.c_str(), &statBuf)==0 && S_ISREG(statBuf.st_mode);
}

void copyIfExists(const string &dirIn, const string &dirOut, const string &fileName)
{
    const string pathIn = dirIn + "/" + fileName;
    const string pathOut = dirOut + "/" + fileName;
    if (pathIn==pathOut) {
        return;
    };

    remove(pathOut.c_str());
    if (fileExists(pathIn)) {
        copyFile(pathIn, pathOut);
    };
}

void copyPreferredIfExists(const string &dirPreferred, const string &dirFallback, const string &dirOut, const string &fileName)
{
    const string pathPreferred = dirPreferred + "/" + fileName;
    const string pathFallback = dirFallback + "/" + fileName;
    const string pathOut = dirOut + "/" + fileName;
    if (fileExists(pathPreferred)) {
        if (pathPreferred!=pathOut) {
            remove(pathOut.c_str());
            copyFile(pathPreferred, pathOut);
        };
    } else if (fileExists(pathFallback)) {
        if (pathFallback!=pathOut) {
            remove(pathOut.c_str());
            copyFile(pathFallback, pathOut);
        };
    } else {
        remove(pathOut.c_str());
    };
}

uint64 parseUint64Token(const string &token, const string &path, Parameters &P)
{
    uint64 value=0;
    istringstream tokenStream(token);
    tokenStream >> value;
    if (tokenStream.fail()) {
        ostringstream errOut;
        errOut << "EXITING because of fatal ERROR: could not parse integer token '" << token << "' while writing " << path << "\n";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_GENOME_FILES, P);
    };
    return value;
}

vector<string> splitFields(const string &line)
{
    vector<string> fields;
    istringstream lineStream(line);
    string field;
    while (lineStream >> field) {
        fields.push_back(field);
    };
    return fields;
}

string joinFields(const vector<string> &fields)
{
    if (fields.size()==0) {
        return "";
    };

    string line=fields.at(0);
    for (uint64 ii=1; ii<fields.size(); ii++) {
        line += "\t" + fields.at(ii);
    };
    return line;
}

void offsetField(vector<string> &fields, const uint64 fieldIndex, const uint64 offset, const string &path, Parameters &P)
{
    if (fields.size()<=fieldIndex) {
        ostringstream errOut;
        errOut << "EXITING because of fatal ERROR: malformed annotation sidecar line while writing " << path << "\n";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_GENOME_FILES, P);
    };
    fields.at(fieldIndex)=to_string(parseUint64Token(fields.at(fieldIndex), path, P)+offset);
}

bool readCountedFile(const string &path, uint64 &count, vector<string> &lines, Parameters &P)
{
    count=0;
    lines.clear();
    if (!fileExists(path)) {
        return false;
    };

    ifstream fileIn(path.c_str());
    if (fileIn.fail()) {
        ostringstream errOut;
        errOut << "EXITING because of fatal ERROR: could not open annotation sidecar file " << path << "\n";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_GENOME_FILES, P);
    };

    string line;
    if (!getline(fileIn, line)) {
        return true;
    };
    count=parseUint64Token(line, path, P);

    while (getline(fileIn, line)) {
        lines.push_back(line);
    };

    return true;
}

bool readPlainFile(const string &path, vector<string> &lines, Parameters &P)
{
    lines.clear();
    if (!fileExists(path)) {
        return false;
    };

    ifstream fileIn(path.c_str());
    if (fileIn.fail()) {
        ostringstream errOut;
        errOut << "EXITING because of fatal ERROR: could not open annotation sidecar file " << path << "\n";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_GENOME_FILES, P);
    };

    string line;
    while (getline(fileIn, line)) {
        lines.push_back(line);
    };

    return true;
}

void mergeCountedFile(const string &dirBase, const string &dirInsert, const string &dirOut, const string &fileName, const function<string(const string &)> &transformInsertLine, Parameters &P)
{
    const string pathBase=dirBase + "/" + fileName;
    const string pathInsert=dirInsert + "/" + fileName;
    const string pathOut=dirOut + "/" + fileName;

    uint64 countBase=0, countInsert=0;
    vector<string> linesBase, linesInsert;
    const bool baseExists=readCountedFile(pathBase, countBase, linesBase, P);
    const bool insertExists=readCountedFile(pathInsert, countInsert, linesInsert, P);

    remove(pathOut.c_str());
    if (!baseExists && !insertExists) {
        return;
    };

    ofstream &fileOut=ofstrOpen(pathOut, ERROR_OUT, P);
    fileOut << countBase+countInsert << "\n";
    for (auto &line : linesBase) {
        fileOut << line << "\n";
    };
    for (auto &line : linesInsert) {
        fileOut << transformInsertLine(line) << "\n";
    };
    fileOut.close();
}

void mergePlainFile(const string &dirBase, const string &dirInsert, const string &dirOut, const string &fileName, const function<string(const string &)> &transformInsertLine, Parameters &P)
{
    const string pathBase=dirBase + "/" + fileName;
    const string pathInsert=dirInsert + "/" + fileName;
    const string pathOut=dirOut + "/" + fileName;

    vector<string> linesBase, linesInsert;
    const bool baseExists=readPlainFile(pathBase, linesBase, P);
    const bool insertExists=readPlainFile(pathInsert, linesInsert, P);

    remove(pathOut.c_str());
    if (!baseExists && !insertExists) {
        return;
    };

    ofstream &fileOut=ofstrOpen(pathOut, ERROR_OUT, P);
    for (auto &line : linesBase) {
        fileOut << line << "\n";
    };
    for (auto &line : linesInsert) {
        fileOut << transformInsertLine(line) << "\n";
    };
    fileOut.close();
}

string offsetTranscriptInfoLine(const string &line, const uint64 geneOffset, const uint64 exonOffset, const bool adjustTrend, uint64 &nextTrend, const string &path, Parameters &P)
{
    vector<string> fields=splitFields(line);
    if (adjustTrend) {
        // transcriptInfo.tab carries a rolling transcript-end field. Preserve
        // that convention so genomeInsert output matches a full rebuild.
        if (fields.size()<=3) {
            ostringstream errOut;
            errOut << "EXITING because of fatal ERROR: malformed transcriptInfo.tab line while writing " << path << "\n";
            exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_GENOME_FILES, P);
        };
        const uint64 transcriptEnd=parseUint64Token(fields.at(2), path, P);
        fields.at(3)=to_string(nextTrend);
        nextTrend=max(nextTrend, transcriptEnd);
    };
    offsetField(fields, 6, exonOffset, path, P);
    offsetField(fields, 7, geneOffset, path, P);
    return joinFields(fields);
}

string offsetExonGeTrInfoLine(const string &line, const uint64 geneOffset, const uint64 transcriptOffset, const string &path, Parameters &P)
{
    vector<string> fields=splitFields(line);
    offsetField(fields, 3, geneOffset, path, P);
    offsetField(fields, 4, transcriptOffset, path, P);
    return joinFields(fields);
}

string offsetSjdbListFromGTFLine(const string &line, const uint64 geneOffset, const string &path, Parameters &P)
{
    vector<string> fields=splitFields(line);
    if (fields.size()<5) {
        ostringstream errOut;
        errOut << "EXITING because of fatal ERROR: malformed sjdbList.fromGTF.out.tab line while writing " << path << "\n";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_GENOME_FILES, P);
    };

    string genesOut;
    size_t geneStart=0;
    while (geneStart<=fields.at(4).size()) {
        size_t geneEnd=fields.at(4).find(',', geneStart);
        string geneToken=fields.at(4).substr(geneStart, geneEnd==string::npos ? string::npos : geneEnd-geneStart);
        if (geneToken=="") {
            break;
        };
        if (genesOut!="") {
            genesOut += ",";
        };
        genesOut += to_string(parseUint64Token(geneToken, path, P)+geneOffset);
        if (geneEnd==string::npos) {
            break;
        };
        geneStart=geneEnd+1;
    };
    fields.at(4)=genesOut;

    return joinFields(fields);
}
}

void genomeInsertCopyAnnotationSidecars(const string &dirIn, const string &dirOut)
{
    for (auto &fileName : annotationSidecarFiles) {
        copyIfExists(dirIn, dirOut, fileName);
    };
}

void genomeInsertMergeAnnotationSidecars(const string &dirBase, const string &dirInsert, const string &dirOut, Parameters &P, bool copySjdbFiles)
{
    if (copySjdbFiles) {
        copyPreferredIfExists(dirInsert, dirBase, dirOut, "sjdbInfo.txt");
        copyPreferredIfExists(dirInsert, dirBase, dirOut, "sjdbList.out.tab");
    };

    uint64 geneOffset=0, transcriptOffset=0, exonOffset=0;
    vector<string> unusedLines, transcriptLinesBase;
    readCountedFile(dirBase+"/geneInfo.tab", geneOffset, unusedLines, P);
    readCountedFile(dirBase+"/transcriptInfo.tab", transcriptOffset, transcriptLinesBase, P);
    readCountedFile(dirBase+"/exonInfo.tab", exonOffset, unusedLines, P);

    bool adjustTranscriptTrend=false;
    uint64 nextTranscriptTrend=0;
    if (transcriptLinesBase.size()>0) {
        vector<string> transcriptFields=splitFields(transcriptLinesBase.back());
        if (transcriptFields.size()<=3) {
            ostringstream errOut;
            errOut << "EXITING because of fatal ERROR: malformed transcriptInfo.tab line while writing " << dirOut+"/transcriptInfo.tab" << "\n";
            exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_GENOME_FILES, P);
        };
        nextTranscriptTrend=max(parseUint64Token(transcriptFields.at(2), dirBase+"/transcriptInfo.tab", P),
                                parseUint64Token(transcriptFields.at(3), dirBase+"/transcriptInfo.tab", P));
        adjustTranscriptTrend=true;
    };

    mergePlainFile(dirBase, dirInsert, dirOut, "sjdbList.fromGTF.out.tab",
                   [&](const string &line) { return offsetSjdbListFromGTFLine(line, geneOffset, dirOut+"/sjdbList.fromGTF.out.tab", P); },
                   P);
    mergeCountedFile(dirBase, dirInsert, dirOut, "geneInfo.tab",
                     [&](const string &line) { return line; },
                     P);
    mergeCountedFile(dirBase, dirInsert, dirOut, "transcriptInfo.tab",
                     [&](const string &line) { return offsetTranscriptInfoLine(line, geneOffset, exonOffset, adjustTranscriptTrend, nextTranscriptTrend, dirOut+"/transcriptInfo.tab", P); },
                     P);
    mergeCountedFile(dirBase, dirInsert, dirOut, "exonInfo.tab",
                     [&](const string &line) { return line; },
                     P);
    mergeCountedFile(dirBase, dirInsert, dirOut, "exonGeTrInfo.tab",
                     [&](const string &line) { return offsetExonGeTrInfoLine(line, geneOffset, transcriptOffset, dirOut+"/exonGeTrInfo.tab", P); },
                     P);
}
