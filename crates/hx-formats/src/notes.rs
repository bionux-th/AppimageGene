//! Local Notes field updates. Preserve unknown fields, attributes and rich text
//! outside the field explicitly changed; never rebuild the entire Notes block.
pub fn field(xml:&[u8],tag:&str)->Option<String> {
    let text=std::str::from_utf8(xml).ok()?;
    let (a,b)=content_range(text,tag)?;
    Some(unescape(&text[a..b]))
}

fn content_range(xml:&str,tag:&str)->Option<(usize,usize)> {
    let prefix=format!("<{tag}");
    let mut from=0;
    while let Some(at)=xml[from..].find(&prefix) {
        let at=from+at;let tail=at+prefix.len();
        if !matches!(xml.as_bytes().get(tail),Some(b'>')|Some(b' ')|Some(b'\t')|Some(b'\n')|Some(b'\r')) {from=tail;continue;}
        let end=tail+xml[tail..].find('>')?;
        if xml[..end].ends_with('/') {return None;}
        let a=end+1;let b=a+xml[a..].find(&format!("</{tag}>"))?;return Some((a,b));
    }
    None
}

pub fn set_field(xml:&[u8],tag:&str,value:&str)->Vec<u8> {
    assert!(tag.bytes().all(|b|b.is_ascii_alphanumeric()||b==b'_'));
    let mut text=String::from_utf8_lossy(xml).into_owned();
    let escaped=escape(value);
    if let Some((a,b))=content_range(&text,tag) {text.replace_range(a..b,&escaped);}
    else {
        // Empty fields in native Notes may be self-closing.
        let empty=format!("<{tag}/>");let spaced=format!("<{tag} />");
        let full=format!("<{tag}>{escaped}</{tag}>");
        if text.contains(&empty){text=text.replacen(&empty,&full,1);}
        else if text.contains(&spaced){text=text.replacen(&spaced,&full,1);}
        else if let Some(at)=text.rfind("</Notes>"){text.insert_str(at,&full);}
    }
    text.into_bytes()
}

pub fn escape(s:&str)->String {s.replace('&',"&amp;").replace('<',"&lt;").replace('>',"&gt;").replace('"',"&quot;").replace('\'',"&apos;")}
fn unescape(s:&str)->String {s.replace("&lt;","<").replace("&gt;",">").replace("&quot;","\"").replace("&apos;","'").replace("&amp;","&")}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn updates_only_selected_field_and_keeps_unknown_xml_and_attributes() {
        let xml=br#"<Notes extension="keep"><Created UTC="01:02:03">2026.1.1</Created><Unknown a="b">&amp;</Unknown><CustomMapLabel>Untitled</CustomMapLabel><Description/><UseCustomMapLabel>1</UseCustomMapLabel></Notes>"#;
        let changed=set_field(xml,"CustomMapLabel","445566 & RNA");
        assert_eq!(field(&changed,"CustomMapLabel").as_deref(),Some("445566 & RNA"));
        assert!(String::from_utf8_lossy(&changed).contains("<Unknown a=\"b\">&amp;</Unknown>"));
        assert!(String::from_utf8_lossy(&changed).contains("<Created UTC=\"01:02:03\">2026.1.1</Created>"));
        let changed=set_field(&changed,"Description","test");assert_eq!(field(&changed,"Description").as_deref(),Some("test"));
    }
}
