#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <shobjidl.h>

#include <iostream>

namespace {
int fail(const wchar_t* operation, HRESULT result)
{
    std::wcerr << operation << L" failed (0x" << std::hex << result << L")\n";
    return 1;
}
}

int wmain(int argc, wchar_t* argv[])
{
    if (argc != 4) {
        std::wcerr << L"Usage: create_shortcut <shortcut> <target> <relative-target>\n";
        return 2;
    }

    const HRESULT initialized = CoInitializeEx(nullptr, COINIT_APARTMENTTHREADED);
    if (FAILED(initialized)) return fail(L"CoInitializeEx", initialized);

    IShellLinkW* link = nullptr;
    HRESULT result = CoCreateInstance(CLSID_ShellLink, nullptr, CLSCTX_INPROC_SERVER,
        IID_IShellLinkW, reinterpret_cast<void**>(&link));
    if (FAILED(result)) { CoUninitialize(); return fail(L"CoCreateInstance", result); }

    result = link->SetPath(argv[2]);
    if (SUCCEEDED(result)) result = link->SetRelativePath(argv[3], 0);
    if (SUCCEEDED(result)) {
        std::wstring workingDirectory(argv[2]);
        const auto separator = workingDirectory.find_last_of(L"\\/");
        if (separator != std::wstring::npos) {
            workingDirectory.resize(separator);
            result = link->SetWorkingDirectory(workingDirectory.c_str());
        }
    }
    if (SUCCEEDED(result)) result = link->SetDescription(L"Navier — Transient Model");

    IPersistFile* file = nullptr;
    if (SUCCEEDED(result)) result = link->QueryInterface(IID_IPersistFile, reinterpret_cast<void**>(&file));
    if (SUCCEEDED(result)) result = file->Save(argv[1], TRUE);

    if (file) file->Release();
    link->Release();
    CoUninitialize();
    return FAILED(result) ? fail(L"Creating shortcut", result) : 0;
}
